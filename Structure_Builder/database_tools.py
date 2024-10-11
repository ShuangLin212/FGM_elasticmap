def get_structures_from_database(db, prototype, subl_model, subl_site_ratios):
    """Returns a list of Structure objects from the db that match the criteria.

    The returned list format supports matching SQS to phases that have multiple solution sublattices
    and the inclusion of higher and lower ordered SQS that match the criteria.

    Parameters
    ----------
    db : tinydb.database.Table
        TinyDB database of the SQS database
    prototype : str
        Prototype symbol as in ATAT sqsdb, e.g. 'GAMMA_L12'.
    subl_model : [[str]]
        List of strings of species names. This sublattice model can be of higher dimension than the SQS.
        Outer dimension should be the same length as subl_site_ratios.
    subl_site_ratios : [[float]]
        Scalar multiple of site ratios of each sublattice. e.g. [1, 2] will match [2, 4] and vice
        versa. Outer dimension should be the same length as subl_model.

    Returns
    -------
    [AbstractSQS]
        Abstract SQSs that match the symmetry and sublattice model.
    """
    def lists_are_multiple(l1, l2):
        """
        Returns True if list a is a multiple of b or vice versa

        Parameters
        ----------
        l1 : [int]
        l2 : [int]

        Returns
        -------
        bool

        """
        # can we compare these two lists?
        if len(l1) != len(l2):
            return False
        for a, b in [(l1, l2), (l2, l1)]:  # check the reverse of the lists too
            # see if all a are perfectly divisible by all b
            if all([(x % y == 0) for x, y in zip(a, b)]):
                # see if all have the same multiple
                if len(set([x/y for x, y in zip(a, b)])) == 1:
                    return True
        return False

    from tinydb import where
    results = db.search((where('prototype') == prototype) &
                        (where('sublattice_site_ratios').test(
                         lambda x: (lists_are_multiple([sum(subl) for subl in x], subl_site_ratios))))
              )
    return results