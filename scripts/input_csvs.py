import uproot
import pandas as pd
import glob
import os

root_dir = os.path.join(os.path.dirname(__file__), "../training_data/run1_data/")
root_files = glob.glob(os.path.join(root_dir, "*.root"))

for root_path in root_files:
    base = os.path.basename(os.path.splitext(root_path)[0])
    csv_name = f"{base}.csv"
    csv_path = os.path.join(root_dir, csv_name)
    
    if os.path.exists(csv_path):
        print(f"Skipping {base}.root - CSV already exists")
        continue
    
    try:
        print(f"Creating CSV for {base}.root  ->  {csv_name}")
        f = uproot.open(root_path)
        tree = f["image2d_image2d_binary_tree"]
        
        # Try old format first
        try:
            runs = tree["_run"].array(library="np")
            subruns = tree["_subrun"].array(library="np")
            events = tree["_event"].array(library="np")
        except:
            # Fall back to new nested format
            runs = tree["image2d_image2d_binary_branch/larcv::EventBase/_run"].array(library="np")
            subruns = tree["image2d_image2d_binary_branch/larcv::EventBase/_subrun"].array(library="np")
            events = tree["image2d_image2d_binary_branch/larcv::EventBase/_event"].array(library="np")
        
        df = pd.DataFrame({
            "run_number": runs,
            "subrun_number": subruns,
            "event_number": events
        })
        df.to_csv(csv_path, index=False)
        print(f"   Saved {csv_name}   ({len(df)} entries)")
        
    except Exception as e:
        print(f"SKIPPING {base}.root - Error: {e}")

print("\nDone!")
