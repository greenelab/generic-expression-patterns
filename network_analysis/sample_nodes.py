import string

import numpy as np

def sample_degree_matched(nodes, degrees, is_generic, num_bins=10):

    """
    Sample nodes from a network, with rough degree matching to is_generic.

    Essentially, this works as follows:
      * Split array into num_bins bins, of approximately equal size
        (that is, all the nodes in each bin have approximately equal degree)
      * For each generic node, sample a non-generic node from the same bin,
        without replacement
        (if there are more generic nodes than non-generic nodes in a bin,
         an error will be thrown at this step)

    Arguments
    ---------
    nodes (np.array): node names
    degrees (np.array): node degrees
    is_generic (np.array): boolean array, nodes are or are not generic
    num_bins (int): number of bins to split node list into

    Returns
    -------
    sampled_nodes (np.array): names of sampled nodes
    sampled_degrees (np.array): degrees of sampled arrays
    bin_ixs (np.array): array of arrays containing indexes in each bin

    NOTE: nodes, degrees, is_generic are assumed to be sorted by ascending
    degree when passed into the function.
    """
    # make sure all arrays have equal shape
    assert nodes.shape == degrees.shape
    assert nodes.shape == is_generic.shape
    # make sure degree array is sorted
    assert np.all(np.diff(degrees) >= 0)

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

        # this should always be the case in our experiments, since there are
        # far fewer generic genes than non-generic genes
        # if it's not we need to have bigger bins or otherwise rethink,
        # since sampling w/o replacement won't work
        assert bin_generic.shape[0] <= bin_non_generic.shape[0]

        # sample as many nodes from the non-generic genes in the bin
        # as there are generic genes in the bin, without replacement
        sample_ixs = np.random.choice(bin_non_generic, size=bin_generic.shape[0])

        # and add to the rest of the bin samples
        sampled_nodes = np.concatenate((sampled_nodes, nodes[sample_ixs]))
        sampled_degrees = np.concatenate((sampled_degrees, degrees[sample_ixs]))

    return (sampled_nodes, sampled_degrees, bin_ixs)


def sort_by_degree(nodes, degrees, is_generic):
    """Sort node info together, by ascending degree."""
    nodes, degrees, is_generic = zip(*sorted(
        zip(nodes, degrees, is_generic),
        key=lambda x: x[1]
    ))
    assert np.all(np.diff(degrees) >= 0)
    return (np.array(nodes),
            np.array(degrees),
            np.array(is_generic))


if __name__ == '__main__':
    np.random.seed(1)
    size=50
    degrees = np.random.geometric(p=0.2, size=size)
    nodes = np.array(list(string.printable)[:size])
    is_generic = np.random.choice([0, 1], size=size, p=[0.8, 0.2])
    # sort nodes and degrees by degree
    nodes, degrees, is_generic = sort_by_degree(nodes, degrees, is_generic)
    print(nodes[is_generic.astype('bool')])
    print(degrees[is_generic.astype('bool')])
    s_nodes, s_degrees, _ = sample_degree_matched(nodes, degrees, is_generic)
    print(s_nodes)
    print(s_degrees)

