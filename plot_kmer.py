#!/usr/bin/env python2.7
"""
This script was part of the jmers development. It makes a simple kmer content plot.
"""

#import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def update(data):
    line.set_ydata(data)
    return line,

def data_gen():
    while True:
        yield np.random.rand(10)

def convert(value):
    v = int(value)
    return np.log(v+1)

if __name__ == '__main__':

    import argparse
    
    parser = argparse.ArgumentParser( description = __doc__ )
    
    parser.add_argument("kmer_content", help="kmer content file")
    
    args = parser.parse_args()
    
    with open(args.kmer_content) as f:
        header = f.readline()
        while header:
            kmers = f.readline()
            kmers = map(convert, kmers.split())
            bases = f.readline()
            bases = map(convert, bases.split())
            fig,ax = plt.subplots()
            line = ax.plot(kmers)
            line2 = ax.plot(bases)
            ax.set_ylim(0, 10)
            plt.show()
            header = f.readline()
        
    
    # fig, ax = plt.subplots()
    # line, = ax.plot(np.random.rand(10))
    # ax.set_ylim(0, 1)
    #
    # ani = animation.FuncAnimation(fig, update, data_gen, interval=100)
    # plt.show()
