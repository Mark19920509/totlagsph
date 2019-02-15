#!/usr/bin/env python
from os import rename, listdir

num = 0
for i in range(1,5):
	for j in range(0,100):
		rename("./drop_{}.{:04d}.png".format(i,j), "./drop_{:04d}.png".format(num))
		num += 1

