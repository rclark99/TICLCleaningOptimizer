import numpy as np

# calculate the metrics from validation results
def get_metrics(uproot_file, id):
    tree = uproot_file[f"jetResponseMetrics" + str(id)]["output"]

    meanResp = tree["meanResp"].array()[0]
    rmsResp  = tree["rmsResp"].array()[0]

    return [meanResp, rmsResp]
    # if minimizing both
    # return [1-meanResp, rmsResp]

# read a csv file, return a matrix
def read_csv(filename):
    matrix = np.genfromtxt(filename, delimiter=",", dtype=float)
    if matrix.ndim == 2:
        return np.genfromtxt(filename, delimiter=",", dtype=float)
    return np.array([matrix])
    
# write a matrix to a csv file
def write_csv(filename, matrix):
    np.savetxt(filename, matrix, fmt='%.18f', delimiter=',')