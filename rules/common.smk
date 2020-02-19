def get_resource(rule,resource):
    ''' Get resource info for rule from config file. Return defaults if not found. '''
    try:
        return config["resources"][rule][resource]
    except KeyError:
        return config["resources"]["default"][resource]
