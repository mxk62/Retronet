SCHEMA = {
    "$schema": "http://json-schema.org/draft-04/schema#",
    "title": "Transform databse entry",
    "description": "A retrosynthetic transform",
    "type": "object",
    "properties": {
        "_id": {
            "type": "number",
            "description": "A unique identifier for a transform",
        },
        "reaction_smarts": {
            "type": "string",
            "description": "A transform specification in SMARTS notation",
        },
        "product_smiles": {
            "type": "array",
            "items": {
                "type": "string"
            },
            "description": "List of SMILES representing possible byproducts"
        },
        "popularity": {
            "type": "number",
            "minimum": 0,
            "description": "Number of existing reactions transform can be applied to"
        }
    },
    "required": ["_id", "reaction_smarts"]
}
