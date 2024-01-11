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
    return [np.add(est_value,np.multiply(N-1,np.subtract(est_value,bias))), err]

def get_super_jackknife_blocks(data, block_size, f=lambda x:x):
    NE = len(data)
    N = []
    data_mean = []
    for e in range(NE):
        N.append(int(len(data[e])/block_size))
        data_mean.append(np.mean(data[e],axis=0))
    block_avgs = []
    for e in range(NE):
        for i in range(N[e]):
            block_av=np.copy(data_mean)
            block_av[e] *= N[e]*block_size
            for j in range(block_size):
                block_av[e] -= data[e][i*block_size+j]
            block_av[e] /= (N[e]-1)*block_size
            block_avgs.append(f(block_av))
    return block_avgs