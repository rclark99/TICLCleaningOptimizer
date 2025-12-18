import patatune as optimizer
import subprocess
from utils import get_metrics, write_csv
import numpy as np
import uproot
import argparse
import os

# parsing argument
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--continuing', type=int, action='store')
parser.add_argument('-d', '--default', action='store_true')
parser.add_argument('-p2', '--phase2', action='store_true')
parser.add_argument('-p', '--num_particles', default=5, type=int, action='store')
parser.add_argument('-i', '--num_iterations', default=1, type=int, action='store')
parser.add_argument('-e', '--num_events', default=100, type=int, action='store')
args = parser.parse_args()

num_iterations = args.num_iterations

optimizer.Logger.setLevel('DEBUG')
optimizer.Randomizer.rng = np.random.default_rng(46)

config = 'reconstructionTICLCleaningByBeta.py'

# bounds on parameters
betaContamMin_lb, betaContamMin_ub = 0.5, 4.0
R0_lb,            R0_ub            = 0.02, 0.30
sigmaT_lb,        sigmaT_ub        = 0.02, 0.20
sigmaDR_lb,       sigmaDR_ub       = 0.02, 0.20
sigmaZ_lb,        sigmaZ_ub        = 5.0,  25.0
tAbsCut_lb,       tAbsCut_ub       = 0.05, 0.50
zAbsCut_lb,       zAbsCut_ub       = 10.0, 60.0
tPower_lb,        tPower_ub        = 0.2,  2.0
drPower_lb,       drPower_ub       = 0.2,  2.0
zPower_lb,        zPower_ub        = 0.5,  3.0
wmin_lb,          wmin_ub          = 1e-4, 1e-2

lb = [
    betaContamMin_lb, R0_lb, sigmaT_lb, sigmaDR_lb, sigmaZ_lb,
    tAbsCut_lb, zAbsCut_lb, tPower_lb, drPower_lb, zPower_lb, wmin_lb
]
ub = [
    betaContamMin_ub, R0_ub, sigmaT_ub, sigmaDR_ub, sigmaZ_ub,
    tAbsCut_ub, zAbsCut_ub, tPower_ub, drPower_ub, zPower_ub, wmin_ub
]

working_dir = 'PSO_TICL_CLEANING_JETMETRICS'

def reco_and_validate(params):
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)

    write_csv(f'{working_dir}/parameters.csv', params)

    validation_result = f'{working_dir}/jet_validation_metrics.root'

    subprocess.run([
        'cmsRun', config,
        'nEvents=' + str(args.num_events),
        f'parametersFile={working_dir}/parameters.csv',
        'outputFile=' + validation_result
    ], check=True)

    print('cmsRun', config,
          'nEvents=' + str(args.num_events),
          f'parametersFile={working_dir}/parameters.csv',
          'outputFile=' + validation_result)

    num_particles = len(params)

    with uproot.open(validation_result) as uproot_file:
        population_fitness = np.array(
            [get_metrics(uproot_file, i) for i in range(num_particles)],
            dtype=float
        )

    return population_fitness

# get default metrics
if args.default:
    defaults = [1.5, 0.10, 0.08, 0.08, 12.5, 0.15, 25.0, 0.5, 0.5, 1.5, 0.001]
    print(f'Len defaults {len(defaults)}')

    default_params = [defaults]
    default_metrics = reco_and_validate(default_params)

    write_csv(f'{working_dir}/default.csv',
              [np.concatenate([default_params[0], default_metrics[0]])])

objective = optimizer.Objective(reco_and_validate, 2)
optimizer.FileManager.working_dir = working_dir
optimizer.FileManager.loading_enabled = False
optimizer.FileManager.saving_enabled = True

print(f"Len ub {len(ub)}, lb {len(lb)}")
pso = optimizer.MOPSO(
    objective=objective,
    lower_bounds=lb,
    upper_bounds=ub,
    num_particles=args.num_particles,
    inertia_weight=0.5,
    cognitive_coefficient=1.5,
    social_coefficient=1.5
)

pso.optimize(num_iterations, max_iterations_without_improvement=10)
