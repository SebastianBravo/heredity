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


print(joint_probability(people, {}, {"Lily"}, {"James"}))