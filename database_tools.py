
"""
Functions for querying structures from a database.
"""

from tinydb import where

def get_structures_from_database(db, prototype, subl_model, subl_site_ratios):
    """
    Returns a list of Structure objects from the db that match the criteria.
    """
    def lists_are_multiple(l1, l2):
        """
        Returns True if list a is a multiple of b or vice versa.
        """
        if len(l1) != len(l2):
            return False
        for a, b in [(l1, l2), (l2, l1)]:
            if all([(x % y == 0) for x, y in zip(a, b)]):
                if len(set([x / y for x, y in zip(a, b)])) == 1:
                    return True
        return False

    results = db.search((where('prototype') == prototype) &
                        (where('sublattice_site_ratios').test(
                         lambda x: (lists_are_multiple([sum(subl) for subl in x], subl_site_ratios))))
              )
    return results
