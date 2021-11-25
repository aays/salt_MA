"""
metabolic_knockouts.py - perform single gene knockouts using
iCre1355 mixotrophic network and return predicted biomass
"""

import time
import csv
import cobra
import argparse
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(
        description='single gene knockout simulations', 
        usage='python metabolic_knockouts.py [options]')

    parser.add_argument('-f', '--model_fname', required=True,
        type=str, help='Path to input sbml model')
    parser.add_argument('--objective', required=True,
        type=str, help='Objective to calc biomass with')
    parser.add_argument('-o', '--out', required=True,
        type=str, help='File to write to')

    args = parser.parse_args()

    return args.model_fname, args.objective, args.out

def knockouts(model_fname, objective, out):
    """Perform single gene knock outs and compute biomass

    Biomass should be one of ['mixo', 'auto', 'hetero']

    Script stops for 0.1s after every iteration since cobrapy
    has some weird memory issues (which is why I'm not just
    using the single gene deletion function - it hangs every
    time after ~1000 iterations) and this helps catch where
    it hangs if it does

    Parameters
    ----------
    model_fname : str
        path to sbml model
    objective : str
        objective to calculate (mixo|auto|hetero)
    out : str
        file to write to


    Returns
    -------
    None
    """
    print('[saltMA] loading network...')
    model = cobra.io.read_sbml_model(model_fname)
    model.objective = f'Biomass_Chlamy_{objective}' # mixo, hetero, auto
    print(f'[saltMA] objective is {model.objective}')
    wt_growth = round(model.optimize().objective_value, 3)

    with open(out, 'w', newline='') as f:
        fieldnames = [
            'objective', 'gene', 'biomass_knockout', 'biomass_original',
            'biomass_difference', 'biomass_ratio', 'growth', 'status']
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        print('[saltMA] performing knockouts...')

        for gene in tqdm(model.genes):
            with model:
                gene.knock_out()
                solution = model.optimize()
                obj_val = round(solution.objective_value, 3)
                growth = False if round(obj_val, 2) == 0 else True
                out_dict = {
                    'objective': objective,
                    'gene': gene.id,
                    'biomass_knockout': obj_val,
                    'biomass_original': wt_growth,
                    'biomass_difference': obj_val - wt_growth,
                    'biomass_ratio': obj_val / wt_growth,
                    'growth': growth,
                    'status': solution.status
                    }
                writer.writerow(out_dict)
                time.sleep(0.1)


def main():
    knockouts(*args())

if __name__ == '__main__':
    main()

        

