from heredity import *

people = {
	'Harry': {'name': 'Harry', 'mother': 'Lily', 'father': 'James', 'trait': None},
	'James': {'name': 'James', 'mother': None, 'father': None, 'trait': True},
	'Lily': {'name': 'Lily', 'mother': None, 'father': None, 'trait': False}
}

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

for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")