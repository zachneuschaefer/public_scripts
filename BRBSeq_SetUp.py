import pandas as pd
import string
import tkinter as tk
from tkinter import filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px

def process_rna_file():
    # Open file dialog to select input CSV
    file_path = filedialog.askopenfilename(title="Select RNA CSV File", filetypes=[("CSV Files", "*.csv")])
    if not file_path:
        return
    
    # Load the CSV file
    try:
        rna_df = pd.read_csv(file_path)
    except Exception as e:
        messagebox.showerror("File Error", f"Error reading CSV file: {str(e)}")
        return
    
    # Check for required columns
    required_columns = ["ID", "RIN", "RNA Concentration (ng/uL)", "RNA Source"]
    missing_columns = [col for col in required_columns if col not in rna_df.columns]
    if missing_columns:
        messagebox.showerror("Column Error", f"Missing required columns: {', '.join(missing_columns)}")
        return
    
    # Extract relevant columns
    filtered_rna_df = rna_df[required_columns].copy()

    # Split 'RNA Source' into 'Source Plate' and 'Source Well'
    filtered_rna_df[['Source Plate', 'Source Well']] = filtered_rna_df["RNA Source"].str.split(":", expand=True)

    # Convert RIN and RNA Concentration to numeric values
    filtered_rna_df["RIN"] = pd.to_numeric(filtered_rna_df["RIN"], errors='coerce')
    filtered_rna_df["RNA Concentration (ng/uL)"] = pd.to_numeric(filtered_rna_df["RNA Concentration (ng/uL)"], errors='coerce')

    # Sort by ascending RIN values
    filtered_rna_df = filtered_rna_df.sort_values(by="RIN", ascending=True).reset_index(drop=True)

    # Get user inputs for RNA normalization
    try:
        rna_input_amount = float(rna_input_entry.get())
        final_volume = float(final_volume_entry.get())
        
        if rna_input_amount <= 0 or final_volume <= 0:
            messagebox.showerror("Input Error", "RNA input amount and final volume must be positive values.")
            return
            
    except ValueError:
        messagebox.showerror("Input Error", "Please enter valid numeric values for RNA input and final volume.")
        return

    # Calculate RNA and water volumes
    filtered_rna_df["RNA Volume (uL)"] = rna_input_amount / filtered_rna_df["RNA Concentration (ng/uL)"]
    filtered_rna_df["Water Volume (uL)"] = final_volume - filtered_rna_df["RNA Volume (uL)"]
    filtered_rna_df.loc[filtered_rna_df["Water Volume (uL)"] < 0, "Water Volume (uL)"] = 0

    # Convert volumes to nanoliters
    filtered_rna_df["RNA Volume (nL)"] = filtered_rna_df["RNA Volume (uL)"] * 1000
    filtered_rna_df["Water Volume (nL)"] = filtered_rna_df["Water Volume (uL)"] * 1000

    # Generate destination wells dynamically (384 wells per plate)
    num_samples = len(filtered_rna_df)
    num_dest_plates = (num_samples // 384) + (1 if num_samples % 384 else 0)

    destination_wells = [f"{row}{col}" for plate in range(num_dest_plates) 
                         for col in range(1, 25) for row in string.ascii_uppercase[:16]][:num_samples]

    destination_plates = [f"384-well Eppendorf PCR Plate {i+1}" for i in range(num_dest_plates) for _ in range(384)][:num_samples]

    filtered_rna_df["Destination Well"] = destination_wells
    filtered_rna_df["Destination Plate"] = destination_plates

    # Get output directory - moved this up before we start saving files
    output_dir = filedialog.askdirectory(title="Select Output Directory")
    if not output_dir:
        return
        
    # Save Updated Sample Information
    sample_info_path = f"{output_dir}/Updated_Sample_Information.csv"
    filtered_rna_df.to_csv(sample_info_path, index=False)

    # Create RNA Worklist
    rna_worklist_path = f"{output_dir}/RNA_Worklist.csv"
    rna_worklist_df = filtered_rna_df[["Source Plate", "Source Well", "Destination Plate", "Destination Well", "RNA Volume (nL)"]].copy()
    rna_worklist_df.rename(columns={
        "RNA Volume (nL)": "Transfer Volume"
    }, inplace=True)
    rna_worklist_df.to_csv(rna_worklist_path, index=False)

    # Create H2O Worklist
    h2o_worklist_path = f"{output_dir}/H2O_Worklist.csv"
    water_sources = ["A01", "A02", "A03", "B01", "B02", "B03"]
    h2o_worklist_df = filtered_rna_df.copy()
    h2o_worklist_df["Source Plate"] = "6RES_AQ_BP2"
    h2o_worklist_df["Source Well"] = [water_sources[i % len(water_sources)] for i in range(len(filtered_rna_df))]
    h2o_worklist_df = h2o_worklist_df[["Source Plate", "Source Well", "Destination Plate", "Destination Well", "Water Volume (nL)"]].copy()
    h2o_worklist_df.rename(columns={
        "Water Volume (nL)": "Transfer Volume"
    }, inplace=True)
    h2o_worklist_df.to_csv(h2o_worklist_path, index=False)

    # First Strand Synthesis Reagents Worklist
    reagents_worklist_path = f"{output_dir}/FirstStrandSynthesisReagents.csv"
    
    reagent_sources = [f"{row}{col:02d}" for row in "ABC" for col in range(1, 25)]
    
    def assign_wells(volume_list, source_wells):
        source_index = 0
        cumulative_vol = 0
        assigned_wells = []
        for vol in volume_list:
            cumulative_vol += vol
            if cumulative_vol > 40000:
                source_index += 1
                cumulative_vol = vol
            assigned_wells.append(source_wells[source_index % len(source_wells)])
        return assigned_wells

    # Replace the existing reagents_worklist_df creation with this updated code
    reagent1_wells = assign_wells([400] * len(filtered_rna_df), reagent_sources[:24])  # Just use A01-A24

    # Explicitly create the reagent2 wells pattern starting with B01
    reagent2_sources = [f"{row}{col:02d}" for row in "BC" for col in range(1, 25)]  # B01-B24, C01-C24
    reagent2_wells = assign_wells([2600] * len(filtered_rna_df), reagent2_sources)

    reagents_worklist_df = pd.DataFrame({
    "Source Plate": ["Reagents_Plate"] * (2 * len(filtered_rna_df)),
    "Source Well": reagent1_wells + reagent2_wells,
    "Destination Plate": filtered_rna_df["Destination Plate"].tolist() * 2,
    "Destination Well": filtered_rna_df["Destination Well"].tolist() * 2,
    "Transfer Volume": [400] * len(filtered_rna_df) + [2600] * len(filtered_rna_df)
})
    reagents_worklist_df.to_csv(reagents_worklist_path, index=False)

    # ERCC Worklist
    ercc_worklist_path = f"{output_dir}/ERCC_Worklist.csv"
    ercc_sources = [f"A{col:02d}" for col in range(1, 25)]
    ercc_wells = assign_wells([100] * len(filtered_rna_df), ercc_sources)

    ercc_worklist_df = pd.DataFrame({
        "Source Plate": ["ERCC"] * len(filtered_rna_df),
        "Source Well": ercc_wells,
        "Destination Plate": filtered_rna_df["Destination Plate"],
        "Destination Well": filtered_rna_df["Destination Well"],
        "Transfer Volume": [100] * len(filtered_rna_df)
    })
    ercc_worklist_df.to_csv(ercc_worklist_path, index=False)

    #Generate Sample Summary Plots
    #Generate Box Plot Distribution of RNA Concentration (ng/uL) for All Samples and Pools
    #Generate Box Plot Distribution of RIN for All Samples and Pools
    rna_df = rna_df.sort_values(by="Library Prep Plate ID")

    #box plot distribution of RIN of All Samples
    fig = px.box(rna_df, y="RIN", points="all")
    fig.show()
    AllRNASamples = fig.write_image(f"{output_dir}/RINofAllSamples.png")

    #box plot distribution of RNA Concentration (ng/uL) of All Samples
    fig = px.box(rna_df, y="RNA Concentration (ng/uL)", points="all")
    fig.show()
    ALLRNAConcentrations = fig.write_image(f"{output_dir}/RNAConcentrationOfAllSamples.png")

    #box plot distributions of RIN across FEX
    fig = px.box(rna_df, x="Fermentation::FEX", y="RIN", color="Fermentation::FEX")
    fig.update_layout(xaxis=dict(title=dict(text='FEX')))
    fig.show()
    ALLRINACROSSFEX = fig.write_image(f"{output_dir}/RINAcrossAllFEX.png")

    #box plot distributions of RNA Concentration (ng/uL) across FEX
    fig = px.box(rna_df, x="Fermentation::FEX", y="RNA Concentration (ng/uL)", color="Fermentation::FEX")
    fig.show()
    ALLRNAConcentrationsACROSSFEX = fig.write_image(f"{output_dir}/RNAConcentrationAcrossAllFEX.png")


    #box plot distributions of RNA Concentration (ng/uL) across Library Prep Plate ID
    fig= px.box(rna_df,x="Library Prep Plate ID", y="RNA Concentration (ng/uL)", points="all", color="Library Prep Plate ID" )
    fig.show()
    RNAConcentrationAcrossPlate = fig.write_image(f"{output_dir}/RNAAcrossPlate.png")

    #box plot distributions of RIN across Library Prep Plate ID
    fig= px.box(rna_df,x="Library Prep Plate ID", y="RIN", points="all", color="Library Prep Plate ID" )
    fig.show()
    RINAcrossPlate = fig.write_image(f"{output_dir}/RINAcrossPlate.png")



    messagebox.showinfo("Success", "All worklists and sample information have been successfully generated!")

# Create GUI window
root = tk.Tk()
root.title("RNA Normalization GUI")

# Labels and entry fields
tk.Label(root, text="Enter RNA Input Amount (ng):").pack(pady=5)
rna_input_entry = tk.Entry(root)
rna_input_entry.pack(pady=5)
rna_input_entry.insert(0, "50")  # Default value

tk.Label(root, text="Enter Final Volume (uL):").pack(pady=5)
final_volume_entry = tk.Entry(root)
final_volume_entry.pack(pady=5)
final_volume_entry.insert(0, "10")  # Default value

# Button to start processing
tk.Button(root, text="Select CSV and Process", command=process_rna_file).pack(pady=20)

# Run GUI
root.mainloop()