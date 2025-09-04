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

def get_jackknife_blocks_2(data1, data2, block_size, f=lambda x:x):
    assert len(data1)=len(data2)
    N = int(len(data1)/block_size)
    data1_mean = np.mean(data1,axis=0)*N*block_size
    data2_mean = np.mean(data2,axis=0)*N*block_size
    block_avgs = []
    for i in range(N):
        block_av1=np.copy(data1_mean)
        block_av2=np.copy(data2_mean)
        for j in range(block_size):
            block_av1 -= data1[i*block_size+j]
            block_av2 -= data2[i*block_size+j]
        block_av1 /= (N-1)*block_size
        block_av2 /= (N-1)*block_size
        block_avgs.append(f(block_av1, block_av2))
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