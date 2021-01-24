from os.path import dirname, join as pjoin
import scipy.io as sio

mat_fname = pjoin('mesh', 'strip2.mat')

print(mat_fname)

mat_contents = sio.loadmat(mat_fname)

sorted(mat_contents.keys())

print(mat_contents['p'].dtype)

print(mat_contents['p'][0][1])

print(mat_contents['p'][1][0])


