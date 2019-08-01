def read_adjusted_coordinates(file):

    """
    This function reads the
    :param file:
    :return:
    """

    stat = ''
    const = ''
    easting = ''
    northing = ''
    zone = ''
    lat = ''
    lon = ''
    h_ortho = ''
    h_ellipse = ''
    x = ''
    y = ''
    z = ''
    sd_e = ''
    sd_n = ''
    sd_u = ''

    lines = []
    go = False
    with open(file) as f:
        for line in f:
            line = line.rstrip()
            if line[:20] == 'Adjusted Coordinates':
                go = True
            if go and line != '':
                lines.append(line)
    for line in lines[5:]:
        print(line)