#NOTE THAT THIS CODE IS ONLY A SLIGHTLY MODIFIED VERSION OF THE ORIGINAL SCRIPT BY HAILEY LOUCKS hloucks@ucsc.edu

import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
import re
import sys

if len(sys.argv) == 1:
   print('No arguments passed')
   sys.exit()

# Function to generate plots from a BED file
def plot_annotations(bed_file):
    # Extract chromosome name from the file name
    chrom_name = re.search(r'chr\d+', bed_file).group(0)
    
    # Load the BED file
    df = pd.read_csv(bed_file, sep="\t", header=None)
    print("READ SUCCESSFULLY")
    
    # Extract relevant columns
    df.columns = ["chrom", "start", "end", "name", "score", "strand", "start2", "end2", "color"]
    print("COLUMNS EXTRACTED")
    
    # Normalize the start and end positions for each chromosome so they start at 0
    df['normalized_start'] = df.groupby('chrom')['start'].transform(lambda x: x - x.min())
    df['normalized_end'] = df['normalized_start'] + (df['end'] - df['start'])
    
    # Split the RGB color codes into lists of integers
    df['color'] = df['color'].apply(lambda x: [int(c) / 255 for c in x.split(',')])
    
    # Determine the maximum normalized end for x-axis limits
    max_normalized_end = df['normalized_end'].max()
    
    # Set x-axis limits
    x_min = 0
    x_max = max_normalized_end
    print(x_max)
    
    # Plotting the annotations
    fig, ax = plt.subplots(figsize=(20, len(df['chrom'].unique()) * 0.5))
    
    chrom_labels = df['chrom'].unique()[::-1]  # Reverse the order of chromosomes
    y_pos_dict = {chrom: i for i, chrom in enumerate(chrom_labels)}  # Create a dictionary for y positions
    
    for chrom in chrom_labels:
        chrom_data = df[df['chrom'] == chrom]
        y_pos = y_pos_dict[chrom]  # Get the correct y position from the dictionary
        for _, row in chrom_data.iterrows():
            start = row['normalized_start']
            end = row['normalized_end']
            color = row['color']
            ax.barh(y_pos, end - start, left=start, color=color)
    
    # Set the y-ticks with chromosome names and increase font size
    ax.set_yticks(range(len(chrom_labels)))
    ax.set_yticklabels(chrom_labels, fontsize=14)
    
    # Adjust the y-axis limits to reduce space at the top and bottom
    ax.set_ylim(-0.5, len(chrom_labels) - 0.5)
    
    # Set plot title with the chromosome name and larger font size
    ax.set_title(f"Annotations for {chrom_name}", fontsize=18)
    
    # Set x-axis limits based on the calculated min and max
    ax.set_xlim(x_min, x_max)
    
    # Remove x-ticks and x-tick labels
    ax.xaxis.set_ticks([])
    ax.xaxis.set_ticklabels([])
    
    # Save the plot without the key
    plt.savefig(f"acrocentrics.{chrom_name}.png", dpi=800, bbox_inches='tight')
    
    # Creating the legend
    color_table = {
        "rDNA": "0,0,0",
        "bSat": "250,153,255",
        "gSat": "172,51,199",
        "CenSat": "0,204,204",
        "activeHOR": "153,0,0",
        "inactiveHOR": "255,102,0",
        "dHOR": "255,146,0",
        "monomeric": "255,204,153",
        "HSat1A": "0,222,96",
        "HSat1B": "27,153,139",
        "HSat2": "0,128,250",
        "HSat3": "51,81,137",
        "ct": "224,224,224",
        "GAP": "169,169,169",
        "mixedAlpha": "204,0,0",
        "HSat2_3": "120,161,187"
    }
    
    legend_elements = [mpatches.Patch(facecolor=[int(x)/255 for x in color.split(',')], label=element)
                       for element, color in color_table.items()]
    
    # Plot the legend
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.legend(handles=legend_elements, loc='center', frameon=False)
    ax.axis('off')
    plt.savefig(f"key.acrocentrics.{chrom_name}.png", dpi=300, bbox_inches='tight')

# Example usage
plot_annotations(sys.argv[1])
