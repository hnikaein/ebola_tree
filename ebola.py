import itertools
import math
import re
from io import StringIO
from multiprocessing.dummy import Pool  # TODO change dummy to pool to use multicore

import numpy as np
from Bio import Phylo
from Bio.Phylo.Consensus import majority_consensus
from ete3 import Tree, TreeStyle, TextFace

from sequence import Sequence
from utils import neighbor_join, upgma

ebola_names = ["Bundibugyo", "Reston", "Sudan", "TaiForest", "Zaire"]
ebola_files = map(lambda x: "ProjectFiles/" + x + "_genome.fasta", ebola_names)
gen_poses = {"NP": (0, 4000), "VP35": (2500, 5000), "VP40": (4000, 6500), "GP": (5500, 9000), "VP30": (8000, 10500),
             "VP24": (10000, 12000), "L": (11000, 19000)}
compute_genome_gen_poses_pre_computed = \
    [('Bundibugyo', 'NP', (581, 418, 3314)), ('Bundibugyo', 'VP35', (246, 3007, 4633)),
     ('Bundibugyo', 'VP40', (185, 4320, 5875)), ('Bundibugyo', 'GP', (426, 5919, 8867)),
     ('Bundibugyo', 'VP30', (187, 8370, 9652)), ('Bundibugyo', 'VP24', (197, 10134, 11504)),
     ('Bundibugyo', 'L', (1792, 11381, 18737)),
     ('Reston', 'NP', (600, 417, 3348)), ('Reston', 'VP35', (248, 3007, 4659)),
     ('Reston', 'VP40', (229, 4350, 5756)), ('Reston', 'GP', (447, 5946, 8859)),
     ('Reston', 'VP30', (194, 8373, 9614)), ('Reston', 'VP24', (247, 10051, 11458)),
     ('Reston', 'L', (1910, 11327, 18794)),
     ('Sudan', 'NP', (558, 444, 3227)), ('Sudan', 'VP35', (292, 3012, 4628)),
     ('Sudan', 'VP40', (232, 4299, 5781)), ('Sudan', 'GP', (427, 5882, 8760)),
     ('Sudan', 'VP30', (209, 8326, 9612)), ('Sudan', 'VP24', (238, 10087, 11474)),
     ('Sudan', 'L', (1828, 11337, 18763)),
     ('Tai Forest', 'NP', (540, 481, 3243)), ('Tai Forest', 'VP35', (237, 3025, 4694)),
     ('Tai Forest', 'VP40', (170, 4303, 5753)), ('Tai Forest', 'GP', (469, 5934, 8733)),
     ('Tai Forest', 'VP30', (222, 8323, 9672)), ('Tai Forest', 'VP24', (226, 10153, 11438)),
     ('Tai Forest', 'L', (1822, 11371, 18861)),
     ('Zaire', 'NP', (562, 449, 3287)), ('Zaire', 'VP35', (243, 3031, 4614)),
     ('Zaire', 'VP40', (270, 4389, 5828)), ('Zaire', 'GP', (424, 5940, 8738)),
     ('Zaire', 'VP30', (192, 8376, 9627)), ('Zaire', 'VP24', (210, 10093, 11495)),
     ('Zaire', 'L', (1908, 11385, 18761))]
compute_genome_gen_distances_pre_computed = \
    [('NP', 0, 1, 1043), ('NP', 0, 2, 1062), ('NP', 0, 3, 880), ('NP', 0, 4, 970), ('NP', 1, 2, 1060),
     ('NP', 1, 3, 1093), ('NP', 1, 4, 1044), ('NP', 2, 3, 1008), ('NP', 2, 4, 991), ('NP', 3, 4, 932),
     ('VP35', 0, 1, 626), ('VP35', 0, 2, 609), ('VP35', 0, 3, 572), ('VP35', 0, 4, 585), ('VP35', 1, 2, 611),
     ('VP35', 1, 3, 667), ('VP35', 1, 4, 617), ('VP35', 2, 3, 676), ('VP35', 2, 4, 603), ('VP35', 3, 4, 598),
     ('VP40', 0, 1, 580), ('VP40', 0, 2, 556), ('VP40', 0, 3, 474), ('VP40', 0, 4, 557), ('VP40', 1, 2, 530),
     ('VP40', 1, 3, 516), ('VP40', 1, 4, 539), ('VP40', 2, 3, 539), ('VP40', 2, 4, 567), ('VP40', 3, 4, 551),
     ('GP', 0, 1, 1096), ('GP', 0, 2, 1152), ('GP', 0, 3, 953), ('GP', 0, 4, 1084), ('GP', 1, 2, 1136),
     ('GP', 1, 3, 1161), ('GP', 1, 4, 1162), ('GP', 2, 3, 1139), ('GP', 2, 4, 1159), ('GP', 3, 4, 974),
     ('VP30', 0, 1, 508), ('VP30', 0, 2, 495), ('VP30', 0, 3, 441), ('VP30', 0, 4, 435), ('VP30', 1, 2, 529),
     ('VP30', 1, 3, 530), ('VP30', 1, 4, 476), ('VP30', 2, 3, 541), ('VP30', 2, 4, 505), ('VP30', 3, 4, 476),
     ('VP24', 0, 1, 566), ('VP24', 0, 2, 549), ('VP24', 0, 3, 468), ('VP24', 0, 4, 513), ('VP24', 1, 2, 554),
     ('VP24', 1, 3, 546), ('VP24', 1, 4, 557), ('VP24', 2, 3, 540), ('VP24', 2, 4, 553), ('VP24', 3, 4, 503),
     ('L', 0, 1, 2461), ('L', 0, 2, 2450), ('L', 0, 3, 1982), ('L', 0, 4, 2197), ('L', 1, 2, 2480), ('L', 1, 3, 2454),
     ('L', 1, 4, 2425), ('L', 2, 3, 2445), ('L', 2, 4, 2399), ('L', 3, 4, 2212)]
compute_genome_distances_pre_computed = [(0, 1, 6797), (0, 2, 6815), (0, 3, 5530), (0, 4, 6247), (0, 5, 9064),
                                         (1, 2, 6788), (1, 3, 6783), (1, 4, 6697), (1, 5, 8926), (2, 3, 6804),
                                         (2, 4, 6740), (2, 5, 8989), (3, 4, 6209), (3, 5, 9046), (4, 5, 8986)]


class MyTreeStyle(TreeStyle):
    """
    A costumized TreeStyle using for fast use in Tree in ete3 package. 
    """

    def __init__(self, title):
        """
        initialize with a title
        """
        super().__init__()
        self.show_branch_length = True
        self.title[1] = [TextFace(title, fsize=20)]


def genome_gen_align(genome, gen):
    """
    align a gen in a genome, and find the best position for that gen in genome
    :return: genome name, gen name and a tuple containing score, start position and end position 
    """
    return genome.get_name(), gen.get_name(), genome.local_global_align(gen, gen_poses[gen.get_name()])


def compute_genome_gen_poses(genomes, gens):
    """
    gets some genomes and some gens and find position of gens in genomes 
    :return: a dictionary that for each genome and gen
    """
    p = Pool()
    # TODO for computing the result remove 0 otherwise the precomputed value will be used
    pool_result = p.starmap(genome_gen_align, list(itertools.product(genomes, gens))[:0])
    p.close()
    if pool_result:  # check if we computed new value or if it should use precomputed value
        print(pool_result)
    else:
        pool_result = compute_genome_gen_poses_pre_computed
    result = {}
    for genome, gen, genomes_gen_pos in pool_result:
        result[(genome, gen)] = genomes_gen_pos
    return result


def init_pool(all_genome_gen_poses1, all_genomes1):
    """
    a dummy function for initializing pools
    """
    global all_genome_gen_poses, all_genomes
    all_genome_gen_poses, all_genomes = all_genome_gen_poses1, all_genomes1


def genome_gen_distance(gen, ii, jj):
    """
    computing edit distance between ii'th genome and jj'th genome base on gene gen.
    :return: the distance
    """
    if ii >= jj:
        return None
    genome_i_poses = all_genome_gen_poses[(all_genomes[ii].get_name(), gen.get_name())]
    genome_j_poses = all_genome_gen_poses[(all_genomes[jj].get_name(), gen.get_name())]
    return gen.get_name(), ii, jj, -1 * all_genomes[ii].sub_genome(genome_i_poses[1], genome_i_poses[2]).global_align(
        all_genomes[jj].sub_genome(genome_j_poses[1], genome_j_poses[2]), match_score=0)


def compute_genome_gen_distances(genomes, gens, genome_gen_poses):
    """
    compute edit distance between all genomes in genomes based on each gen in gens, and by help of genome_gen_poses    
    :return: a dictionary contains a distance matrix for each gene
    """
    p = Pool(initializer=init_pool, initargs=(genome_gen_poses, genomes))
    # TODO for computing the result remove 0 otherwise the precomputed value will be used
    gen_edit_dists = p.starmap(genome_gen_distance,
                               list(itertools.product(gens, range(len(genomes)), range(len(genomes))))[:0])
    p.close()
    if gen_edit_dists:  # check if we computed new value or if it should use precomputed value
        print(gen_edit_dists)
    else:
        gen_edit_dists = compute_genome_gen_distances_pre_computed

    result = {}  # convert list like result to dict of matrixes
    for gen_edit_dist in gen_edit_dists:
        if gen_edit_dist is None:
            continue
        if gen_edit_dist[0] not in result:
            result[gen_edit_dist[0]] = np.full([len(genomes), len(genomes)], 0)
        result[gen_edit_dist[0]][gen_edit_dist[1]][gen_edit_dist[2]] = gen_edit_dist[3]
        result[gen_edit_dist[0]][gen_edit_dist[2]][gen_edit_dist[1]] = gen_edit_dist[3]
    return result


def genome_distance(ii, jj):
    """
    computing edit distance between ii'th genome and jj'th genome
    :return: the distance
    """
    if ii >= jj:
        return None
    return ii, jj, -1 * all_genomes[ii].global_align(all_genomes[jj], match_score=0)


def compute_genome_distances(genomes):
    """
    compute edit distance between all genomes in genomes    
    :return: a distance matrix
    """
    p = Pool(initializer=init_pool, initargs=(None, genomes))
    genome_edit_dists = p.starmap(genome_distance,
                                  list(itertools.product(range(len(genomes)), range(len(genomes))))[:0])  # TODO
    p.close()
    if genome_edit_dists:  # check if we computed new value or if it should use precomputed value
        print(genome_edit_dists)
    else:
        genome_edit_dists = compute_genome_distances_pre_computed

    result = np.full([len(genomes), len(genomes)], 0)
    for genome_edit_dist in genome_edit_dists:  # convert list like result to a matrix
        if genome_edit_dist is None:
            continue
        result[genome_edit_dist[0]][genome_edit_dist[1]] = genome_edit_dist[2]
        result[genome_edit_dist[1]][genome_edit_dist[0]] = genome_edit_dist[2]
    return result


if __name__ == '__main__':
    # reading genes and genomes
    all_gens = Sequence.read_sequence("ProjectFiles/Marburg_Genes.fasta")
    all_genomes = [Sequence("")] * 0
    all_genomes = list(map(lambda x: Sequence.read_sequence(x)[0], ebola_files))

    # computing exact gen poses in genomes
    all_genome_gen_poses = compute_genome_gen_poses(all_genomes, all_gens)
    # print(all_genome_gen_poses)

    # computing edit distance matrix for each gene
    gen_edit_dist_matrixes = compute_genome_gen_distances(all_genomes, all_gens, all_genome_gen_poses)
    # for key, value in gen_edit_dist_matrixes.items(): print(key, "\n", value)

    # drawing trees and writing files for each genes
    upgma_results = []
    nj_results = []
    for gen_name, gen_edit_dist_matrix in gen_edit_dist_matrixes.items():
        # writing matrix to csv file
        file = open(gen_name + ".csv", "w")
        matrix_str = "\n".join(",".join(map(str, row)) for row in gen_edit_dist_matrix)
        file.write(matrix_str)
        file.close()

        # dist_matrix2 = [gen_edit_dist_matrix[i][:i + 1] for i in range(0, len(ebola_names))]

        # draw tree for each algorithm for each gene and store it in a list for consensus tree
        for trees, func, alg in [(upgma_results, upgma, "UPGMA"), (nj_results, neighbor_join, "NJ")]:
            trees.append(func(ebola_names, gen_edit_dist_matrix))
            Tree(trees[-1]).render("%s_%s.png" % (alg, gen_name),
                                   tree_style=MyTreeStyle("Gene %s, %s" % (gen_name, alg)), h=500)

            # print("Gene %s, %s" % (gen_name, alg))
            # Phylo.draw_ascii(Phylo.read(StringIO(trees[-1]), "newick"))
            # Phylo.draw(DistanceTreeConstructor().upgma(DistanceMatrix(ebola_names, dist_matrix2)))

    # add marburg to genomes to have it for edit distance computing
    all_genomes.extend(Sequence.read_sequence("ProjectFiles/Marburg_genome.fasta"))
    # computing edit distance between all genomes
    genome_distances = compute_genome_distances(all_genomes)
    # print(genome_distances)

    # drawing consensus tree and tree of genome edit distance for each algorithm
    for trees, func, alg in [(upgma_results, upgma, "UPGMA"), (nj_results, neighbor_join, "NJ")]:
        temp_io = StringIO()
        # finding consensus tree
        Phylo.write(majority_consensus([Phylo.read(StringIO(tree), "newick") for tree in trees], 2.9 / 7), temp_io,
                    "newick")
        temp_io.seek(0)
        consensus_tree = re.sub(":[\d.]*?:", ":", temp_io.read())
        Tree(consensus_tree).render("%s_consensus.png" % alg, tree_style=MyTreeStyle(alg + " Consensus"), h=500)
        # print("%s_consensus" % alg)
        # Phylo.draw_ascii(Phylo.read(StringIO(consensus_tree), "newick"))
        genome_tree = func(ebola_names, genome_distances)
        Tree(genome_tree).render("%s.png" % alg, tree_style=MyTreeStyle(alg), h=500)
        # print("%s" % alg)
        # Phylo.draw_ascii(Phylo.read(StringIO(genome_tree), "newick"))

    # drawing tree of all genomes (including marburg) by nj algorithm
    all_genomes_tree_with_marburg = neighbor_join(ebola_names + ["Marburg"], genome_distances)
    Tree(all_genomes_tree_with_marburg).render("NJ+.png", tree_style=MyTreeStyle("NJ with Marburg"), h=500)
    # Phylo.draw_ascii(Phylo.read(StringIO(all_genomes_tree_with_marburg), "newick"))

    # computing years between each pair of genomes
    age = np.full([len(all_genomes), len(all_genomes)], 0)
    for i in range(len(genome_distances)):
        for j in range(len(genome_distances)):
            age[i][j] = -.75 * math.log(1 - 4 / 3 * genome_distances[i][j] / all_genomes[i].size) / 1.9 * 1000
    print(age)
