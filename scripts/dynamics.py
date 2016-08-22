#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import gzip
import matplotlib.pyplot as plt
import input

data_size = input.MAX_Z
energy_data_size = 100000

prefix = input.output_directory_name

fig = plt.figure(figsize=(16, 12))

fd = gzip.open('{0}/conc0.gz'.format(prefix), 'rb')
a = np.frombuffer(fd.read(), dtype=np.float32)
b = np.reshape(a, (-1, data_size))

ax = fig.add_subplot(211)
cax = ax.matshow(b, aspect='auto')
fig.colorbar(cax)

fd = gzip.open('{0}/conc1.gz'.format(prefix), 'rb')
a = np.frombuffer(fd.read(), dtype=np.float32)
b = np.reshape(a, (-1, data_size))

ax = fig.add_subplot(212)
cax = ax.matshow(b, aspect='auto')
fig.colorbar(cax)

'''
fd = gzip.open('{0}/conc2.gz'.format(prefix), 'rb')
a = np.frombuffer(fd.read(), dtype=np.float32)
b = np.reshape(a, (-1, data_size))

ax = fig.add_subplot(325)
cax = ax.matshow(b, aspect='auto')
fig.colorbar(cax)
'''

'''
energy_lim = 600.;

fd = gzip.open('{3}/{0}_{1}_{2}/output/energy0.gz'.format(i1, i2, i3, prefix), 'rb')
a = np.frombuffer(fd.read(), dtype=np.int32)
b = np.reshape(a, (-1, energy_data_size))
b1 = b + 1
b_log = np.log(b1)

ax = fig.add_subplot(322)
ax.set_xlim(0, energy_lim/3000.*energy_data_size)
cax = ax.matshow(b_log, aspect='auto', vmin=0, vmax=5)
fig.colorbar(cax)

fd = gzip.open('{3}/{0}_{1}_{2}/output/energy1.gz'.format(i1, i2, i3, prefix), 'rb')
a = np.frombuffer(fd.read(), dtype=np.int32)
b = np.reshape(a, (-1, energy_data_size))
b1 = b + 1
b_log = np.log(b1)

ax = fig.add_subplot(324)
ax.set_xlim(0, energy_lim/3000.*energy_data_size*56./24.)
cax = ax.matshow(b_log, aspect='auto', vmin=0, vmax=3)
fig.colorbar(cax)

fd = gzip.open('{3}/{0}_{1}_{2}/output/energy2.gz'.format(i1, i2, i3, prefix), 'rb')
a = np.frombuffer(fd.read(), dtype=np.int32)
b = np.reshape(a, (-1, energy_data_size))
b1 = b + 1
b_log = np.log(b1)

ax = fig.add_subplot(326)
ax.set_xlim(0, energy_lim/3000.*energy_data_size*197./69.)
cax = ax.matshow(b_log, aspect='auto', vmin=0, vmax=5)
fig.colorbar(cax)

fd = gzip.open('data/0_0_0/output/ez.gz', 'rb')
a = np.frombuffer(fd.read(), dtype=np.double)
b = np.reshape(a, (-1, data_size))

ax = fig.add_subplot(122)
cax = ax.matshow(b, aspect='auto')
fig.colorbar(cax)
'''

plt.show()
