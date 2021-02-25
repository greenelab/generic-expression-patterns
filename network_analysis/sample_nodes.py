import string

import numpy as np

def sample_degree_matched(nodes,
                          degrees,
                          is_generic,
                          sort=False,
                          num_bins=10):

    # NOTE: nodes, degrees, is_generic are assumed to be sorted
    # by ascending degree

    assert nodes.shape[0] == degrees.shape[0]
    assert nodes.shape[0] == is_generic.shape[0]

    # get split indexes
    bin_ixs = np.array_split(np.arange(nodes.shape[0]), num_bins)

    # for each generic node, sample another node from the same degree bin,
    # without replacement
    sampled_nodes, sampled_degrees = np.array([]).astype('int'), np.array([]).astype('int')
    for bin_ix in bin_ixs:
        # indices in the bin that are generic genes
        bin_generic = bin_ix[is_generic[bin_ix].astype('bool')]
        # indices in the bin that are not generic genes
        bin_non_generic = bin_ix[~is_generic[bin_ix].astype('bool')]

        # this should always be the case since there are far fewer
        # generic genes than non-generic genes
        # if it's not we need to have bigger bins or otherwise rethink,
        # since sampling w/o replacement won't work
        assert bin_generic.shape[0] < bin_non_generic.shape[0]

        # sample as many nodes from the non-generic genes in the bin
        # as there are generic genes in the bin, without replacement
        sample_ixs = np.random.choice(bin_non_generic, size=bin_generic.shape[0])

        # and add to the rest of the bin samples
        sampled_nodes = np.concatenate((sampled_nodes, nodes[sample_ixs]))
        sampled_degrees = np.concatenate((sampled_degrees, degrees[sample_ixs]))

    return (sampled_nodes, sampled_degrees, bin_ixs)

if __name__ == '__main__':
    np.random.seed(1)
    size=50
    degrees = np.random.geometric(p=0.2, size=size)
    nodes = np.array(list(string.printable)[:size])
    # sort nodes and degrees by degree
    print(nodes)
    print(degrees)
    nodes, degrees = zip(*sorted(
        zip(nodes, degrees),
        key=lambda x: x[1]
    ))
    nodes = np.array(nodes)
    degrees = np.array(degrees)
    print(nodes)
    print(degrees)
    is_generic = np.random.choice([0, 1], size=size, p=[0.8, 0.2])
    print(is_generic)
    s_nodes, s_degrees, _ = sample_degree_matched(nodes, degrees, is_generic, sort=True)
    print(nodes[is_generic.astype('bool')])
    print(degrees[is_generic.astype('bool')])
    print(s_nodes)
    print(s_degrees)

