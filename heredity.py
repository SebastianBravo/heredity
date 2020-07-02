import csv
import itertools
import sys
import numpy

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    # Initialize an empty set for people without the gene
    no_gene = set()

    # People without parents set
    no_parents_people = no_parents(people)

    # Initialize an empty list to keep track of the probabilities to multiply
    probabilities = []

    # Add people with zero genes to the no_gene set
    for person in people.keys():
        if person not in one_gene and person not in two_genes:
            no_gene.add(person)

    # Loop through the people with zero genes and add their probability to the list
    for person in no_gene:
        probability = 1

        if person in have_trait:
            probability = PROBS["trait"][0][True]
        else:
            probability = PROBS["trait"][0][False]

        probabilities.append(probability)

        if person in no_parents_people:
            probability = PROBS["gene"][0]
        else:
            mother = people[person]["mother"]
            father = people[person]["father"]
            if mother in no_gene:
                if father in no_gene:
                    probability = (1 - PROBS["mutation"]) * (1 - PROBS["mutation"])
                elif father in one_gene:
                    probability = (1 - PROBS["mutation"]) * 0.5
                else:
                    probability = (1 - PROBS["mutation"]) * PROBS["mutation"]
            elif mother in one_gene:
                if father in no_gene:
                    probability = 0.5 * (1 - PROBS["mutation"])
                elif father in one_gene:
                    probability = 0.5 * 0.5
                else:
                    probability = 0.5 * PROBS["mutation"]
            elif mother in two_genes:
                if father in no_gene:
                    probability = PROBS["mutation"] * (1 - PROBS["mutation"])
                elif father in one_gene:
                    probability = PROBS["mutation"] * 0.5
                else:
                    probability = PROBS["mutation"] * PROBS["mutation"]

        probabilities.append(probability)

    # Loop through the people with one gene and add their probability to the list
    for person in one_gene:
        probability = 1

        if person in have_trait:
            probability = PROBS["trait"][1][True]
        else:
            probability = PROBS["trait"][1][False]

        probabilities.append(probability)

        if person in no_parents_people:
            probability = PROBS["gene"][1]
        else:
            mother = people[person]["mother"]
            father = people[person]["father"]
            if mother in no_gene:
                if father in no_gene:
                    probability = 2 * PROBS["mutation"] * (1 - PROBS["mutation"])
                elif father in one_gene:
                    probability = 0.5
                else:
                    probability = 0.9802
            elif mother in one_gene:
                if father in no_gene:
                    probability = 0.5
                elif father in one_gene:
                    probability = 0.5
                else:
                    probability = 0.5
            elif mother in two_genes:
                if father in no_gene:
                    probability = 0.9802
                elif father in one_gene:
                    probability = 0.5
                else:
                    probability = 2 * PROBS["mutation"] * (1 - PROBS["mutation"])

        probabilities.append(probability)

    # Loop through the people with two genes and add their probability to the list
    for person in two_genes:
        probability = 1

        if person in have_trait:
            probability = PROBS["trait"][2][True]
        else:
            probability = PROBS["trait"][2][False]

        probabilities.append(probability)

        if person in no_parents_people:
            probability = PROBS["gene"][2]
        else:
            mother = people[person]["mother"]
            father = people[person]["father"]
            if mother in no_gene:
                if father in no_gene:
                    probability = PROBS["mutation"] * PROBS["mutation"]
                elif father in one_gene:
                    probability = PROBS["mutation"] * 0.5
                else:
                    probability = PROBS["mutation"] * (1 - PROBS["mutation"])
            elif mother in one_gene:
                if father in no_gene:
                    probability = 0.5 * PROBS["mutation"]
                elif father in one_gene:
                    probability = 0.5 * 0.5
                else:
                    probability = 0.5 * (1 - PROBS["mutation"])
            elif mother in two_genes:
                if father in no_gene:
                    probability = (1 - PROBS["mutation"]) * PROBS["mutation"]
                elif father in one_gene:
                    probability = (1 - PROBS["mutation"]) * 0.5
                else:
                    probability = (1 - PROBS["mutation"]) * (1 - PROBS["mutation"])

        probabilities.append(probability)

    # Return the product of all probabilities
    return numpy.product(list(probabilities))


def no_parents(people):
    """
    Return a set with the people that have no parents
    """
    no_parents_people = set()

    for person in people.keys():
        if people[person]["mother"] == None:
            no_parents_people.add(person)

    return no_parents_people

def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    # Loop through all the people and update the probabilities
    for person in probabilities.keys():
        if person in one_gene:
            probabilities[person]["gene"][1] += p
        elif person in two_genes:
            probabilities[person]["gene"][2] += p
        else:
            probabilities[person]["gene"][0] += p

        if person in have_trait:
            probabilities[person]["trait"][True] += p
        else:
            probabilities[person]["trait"][False] += p 


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    for person in probabilities.keys():
        sum_gene = sum(probabilities[person]["gene"].values())
        sum_trait = sum(probabilities[person]["trait"].values())
        alpha_gene = 1 / sum_gene
        alpha_trait = 1 / sum_trait

        for gene in probabilities[person]["gene"].keys():
            probabilities[person]["gene"][gene] *= alpha_gene

        for trait in probabilities[person]["trait"].keys():
            probabilities[person]["trait"][trait] *= alpha_trait


if __name__ == "__main__":
    main()
