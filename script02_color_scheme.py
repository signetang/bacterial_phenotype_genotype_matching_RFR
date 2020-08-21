import matplotlib
matplotlib.use('Agg')
import pandas as pd
import seaborn as sns
sns.set()

def save_color_scheme(color_palette, taxonomy_file, color_file):
    """Use the selected color scheme to make a dictionary of (sub)species to color,
    and save as a text file to be read and used in other scripts."""

    # Read taxonomy data:
    taxonomy = pd.read_csv(taxonomy_file, sep='; ', header=None, usecols=[10], engine='python')
    taxonomy.columns = ['Subspecies']
    # Unknown subspecies:
    taxonomy = taxonomy.fillna('Subspecies unknown')
    # If there's no value in the last column, then there is only a ';' after the last value not a '; ' (which is what is used as sep in pd.read_csv):
    taxonomy['Subspecies'] = taxonomy['Subspecies'].map(lambda x: x.rstrip(';'))

    items = set(taxonomy['Subspecies'])
    n = len(list(items))

    # Use the given list of colors or the color palette:
    if len(color_palette) < 2:
        sns_color_palette = sns.color_palette(color_palette[0], n)
    else:
        sns_color_palette = sns.color_palette(color_palette)
    i = 0
    f = open(color_file,'w')
    for item in items:
        f.write(item)
        for c in sns_color_palette[i]:
            f.write('\t'+str(c))
        f.write('\n')
        i += 1

    f.close()

def load_color_scheme(color_file):
    """Loads the color assignments from the color_file.
    Returns a dictionary of (sub)species to color.
    This function is not run if the script is run as main."""
    colors = {}
    with open(color_file,'r') as f:
        for line in f:
            key,r,b,g = line.split('\t')
            r = float(r)
            b = float(b)
            g = float(g.rstrip('\n'))
            colors[key] = (r,b,g)
    return colors

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--color_palette', help='Seaborn color palette name (for instance "hls") or color names in a list (for instance --color_palette "green,blue,magenta"). No spaces allowed after the commas in the list.')
    parser.add_argument('--taxonomy_file', help="Name and path of the file 'selected_files_based_on_KmerFinder_and_Quast.taxonomy'")
    parser.add_argument('--color_file', help='Path and name of file to save the color assignments to.')
    args = parser.parse_args()

    save_color_scheme(args.color_palette.split(','), args.taxonomy_file, args.color_file)

