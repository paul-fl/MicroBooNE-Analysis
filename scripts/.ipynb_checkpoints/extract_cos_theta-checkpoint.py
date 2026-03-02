import uproot
import pandas as pd
import numpy as np
from pathlib import Path

# NuMI beam direction
v_numi = np.array([0.462372, 0.0488541, 0.885339])

def extract_cos_theta(filepath, tree_path, output_csv):
    """Extract cos(theta) and event info from ntuple."""
    print(f"Processing {filepath}...")
    f = uproot.open(filepath)
    t = f[tree_path]
    
    # Get arrays
    run = t['run_number'].array()
    subrun = t['subrun_number'].array()
    event = t['event_number'].array()
    dirx = t['reco_shower_dirx'].array()
    diry = t['reco_shower_diry'].array()
    dirz = t['reco_shower_dirz'].array()
    
    # Get weights if available
    try:
        weight_spline = t['genie_spline_weight'].array()
        weight_tune = t['genie_CV_tune_weight'].array()
        has_weights = True
    except:
        has_weights = False
    
    # Calculate cos(theta) for each event
    data = []
    for i in range(len(run)):
        if len(dirx[i]) > 0:
            cos_theta = dirx[i][0]*v_numi[0] + diry[i][0]*v_numi[1] + dirz[i][0]*v_numi[2]
            row = {
                'run': int(run[i]), 
                'subrun': int(subrun[i]), 
                'event': int(event[i]), 
                'cos_theta': cos_theta
            }
            if has_weights:
                row['weight_spline'] = weight_spline[i]
                row['weight_tune'] = weight_tune[i]
                row['weight'] = weight_spline[i] * weight_tune[i]
            data.append(row)
    
    df = pd.DataFrame(data)
    df.to_csv(output_csv, index=False)
    print(f"  Saved {len(df)} events to {output_csv}")
    return df

# Output directory
out_dir = Path("/home/paul/Msci/ntuples/cos_theta_csvs/")
out_dir.mkdir(exist_ok=True)

# Signal files
signal_files = list(Path("/home/paul/Msci/ntuples/").glob("run1_dt_ratio_0.6*.root"))
for f in signal_files:
    name = f.stem
    extract_cos_theta(f, 'singlephotonana/vertex_tree', out_dir / f"{name}_cos_theta.csv")

# Nu overlay
extract_cos_theta("/home/paul/Msci/ntuples/run1_NuMI_nu_overlay_with_weights.root", 
                  'vertex_tree', 
                  out_dir / "run1_nu_overlay_cos_theta.csv")

# Dirt
extract_cos_theta("/home/paul/Msci/ntuples/run1_NuMI_dirt_with_weights.root", 
                  'vertex_tree', 
                  out_dir / "run1_dirt_cos_theta.csv")

# Offbeam
extract_cos_theta("/home/paul/Msci/ntuples/run1_NuMI_offbeam_full_set_sp.root", 
                  'singlephotonana/vertex_tree', 
                  out_dir / "run1_offbeam_cos_theta.csv")

# NCPi0
extract_cos_theta("/home/paul/Msci/ntuples/NCPi0_overlay_production_v2_sp.root", 
                  'singlephotonana/vertex_tree', 
                  out_dir / "NCPi0_cos_theta.csv")

# NCeta
extract_cos_theta("/home/paul/Msci/ntuples/NCeta_overlay_production_v2_sp.root", 
                  'singlephotonana/vertex_tree', 
                  out_dir / "NCeta_cos_theta.csv")

print("\nDone!")
