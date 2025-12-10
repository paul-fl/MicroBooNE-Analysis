import uproot
import pandas as pd
import glob
import os

# path to the folder where your ROOT files are
root_dir = os.path.join(os.path.dirname(__file__), "../training_data/run3/")

# find ALL .root files in that directory
root_files = glob.glob(os.path.join(root_dir, "*.root"))
# root_files = [os.path.join(root_dir, "run1_offbeam_larcv_cropped_full_set.root")]


for root_path in root_files:
    base = os.path.basename(os.path.splitext(root_path)[0])
    csv_name = f"{base}.csv"
    csv_path = os.path.join(root_dir, csv_name)

    print(f"Creating CSV for {base}.root  ->  {csv_name}")

    f = uproot.open(root_path)
    tree = f["image2d_image2d_binary_tree"]

    runs     = tree["_run"].array(library="np")
    subruns  = tree["_subrun"].array(library="np")
    events   = tree["_event"].array(library="np")

    df = pd.DataFrame({
        "run_number": runs,
        "subrun_number": subruns,
        "event_number": events
    })

    df.to_csv(csv_path, index=False)
    print(f"   Saved {csv_name}   ({len(df)} entries)")

print("\nAll CSV files created successfully!\n")

