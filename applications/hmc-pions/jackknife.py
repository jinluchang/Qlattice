import numpy as np

def get_jackknife_blocks(data, block_size, f=lambda x:x):
    N = int(len(data)/block_size)
    data_mean = np.mean(data,axis=0)*N*block_size
    block_avgs = []
    for i in range(N):
        block_av=np.copy(data_mean)
        for j in range(block_size):
            block_av -= data[i*block_size+j]
        block_av /= (N-1)*block_size
        block_avgs.append(f(block_av))
    return block_avgs
    
def get_errors_from_blocks(est_value, blocks):
    N = len(blocks)
    err = 0
    bias = 0
    for i in range(N):
        err = np.add(err, (N-1)/N*np.power(np.subtract(est_value,blocks[i]),2))
        bias = np.add(bias,np.divide(blocks[i],N))
    err = np.power(err,0.5)
    bias_correction = np.multiply(N-1,np.subtract(est_value,bias))
    if(np.any(err==0)):
        print(f"Warning: error zero")
    elif(np.any(abs(bias_correction / err) > 1.0)):
        print(f"Warning: Large bias correction {bias_correction} to central value {est_value} with error {err}.")
    return [np.add(est_value,bias_correction), err]