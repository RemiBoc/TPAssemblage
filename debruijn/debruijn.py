#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""
from random import randint
import random
import statistics

import argparse
import os
import sys
import networkx as nx
import matplotlib.pyplot as plt
#from operator import itemgetter
random.seed(9001)

__author__ = "Remi Bocquet"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Remi Bocquet"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Remi Bocquet"
__email__ = "bocquetrem@eisti.eu"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Reads a fastq file and stores it in a generator.
    """
    isfile(fastq_file)
    with open(fastq_file) as fast_q:
        for line in fast_q:
            yield next(fast_q).rstrip("\n")
            next(fast_q)
            next(fast_q)


def cut_kmer(read, kmer_size):
    """Sequences read string with 1-offset between each sequence.
    """
    for i in range(len(read) - kmer_size + 1):
        yield read[i:(i + kmer_size)]


def build_kmer_dict(fastq_file, kmer_size):
    """Builds a dictionary of each sequence created with cut_kmer and each sequence count.
    """
    kmer_dict = {}
    for sequence in read_fastq(fastq_file):
        for kmer in cut_kmer(sequence, kmer_size):
            try:
                kmer_dict[kmer] += 1
            except KeyError:
                kmer_dict[kmer] = 1
    return kmer_dict



def build_graph(kmer_dict):
    """Builds DeBruijn graph from a kmer_dict.
    """
    graph = nx.DiGraph()
    for key in kmer_dict:
        graph.add_edge(key[:-1], key[1:], weight=kmer_dict[key])

    nx.draw(graph, with_labels=True, font_weight='bold')
    plt.show()
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Removes pathes in path_list from graph, and removes starting and sink nodes if required.
    """
    for path in path_list:
        for i in range(1, len(path) - 1):
            graph.remove_node(path[i])

    if delete_entry_node:
        for path in path_list:
            graph.remove_node(path[0])
    if delete_sink_node:
        for path in path_list:
            graph.remove_node(path[-1])

    return graph

def std(data):
    """Returns standard deviation of data.
    """
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """Selects the best path regarding total weight, then length and removes other paths from graph.
    If two path are equivalent, the best path is selected randomly
    """
    selected_paths = []
    selected_paths_lengths = []

    max_weight = max([path_length[i] * weight_avg_list[i] for i in range(len(path_length))])

    for i in range(len(path_list)):
        if path_length[i] * weight_avg_list[i] == max_weight:
            selected_paths.append(path_list[i])
            selected_paths_lengths.append(path_length[i])

    if len(selected_paths) > 1:
        min_len = min(selected_paths_lengths)
        selected_paths = [selected_paths[i] for i in range(len(selected_paths)) if selected_paths_lengths[i] == min_len]

        if (len(selected_paths)) > 1:
            random_int = random.randint(len(selected_paths)) - 1

            path_list.remove(selected_paths[random_int])

    remove_paths(graph, path_list, delete_entry_node, delete_sink_node)

    return graph


def path_average_weight(graph, path):
    """Returns average value of weights of edges in a path.
    """
    return sum([graph[path[i]][path[i+1]]['weight'] for i in range(len(path)-1)]) / (len(path) - 1)


def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    """Returns graph starting nodes.
    """
    return [node for node, in_degree in graph.in_degree() if in_degree == 0]


def get_sink_nodes(graph):
    """Returns graph sink nodes.
    """
    return [node for node, out_degree in graph.out_degree() if out_degree == 0]

def get_contigs(graph, starting_nodes, ending_nodes):
    """Returns graph contigs.
    """
    contigs_list = []
    for start in starting_nodes:
        for end in ending_nodes:
            try:
                path = nx.shortest_path(graph, source=start, target=end)
                contig = path[0]
                for kmer in range(1, len(path)):
                    contig += path[kmer][-1]
                contigs_list.append((contig, len(contig)))
            except nx.NetworkXNoPath:
                continue
    return contigs_list


def save_contigs(contigs_list, output_file):
    """Saves contigs into output_file.
    """
    with open(output_file, 'w') as out_file:
        contig_number = 0
        for contig, contig_length in contigs_list:
            contig_number += 1
            out_file.write(">contig_" + str(contig_number) + " len=" + str(contig_length) + "\n")
            out_file.write(fill(contig) + "\n")


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    data = read_fastq(args.fastq_file)

    graph = build_graph(build_kmer_dict(args.fastq_file, args.kmer_size))

    # Non processed contigs

    contigs_list = get_contigs(graph, get_starting_nodes(graph), get_sink_nodes(graph))
    print(contigs_list)


if __name__ == '__main__':
    main()
