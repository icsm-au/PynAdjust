# pynadjust module for classes used across adj, apu and xyz files

import geodepy.convert as gc


class Station(object):
    def __init__(self, name, description, con, lat, lon, ohgt, ehgt, sd_e, sd_n, sd_u, hpu, vpu, smaj, smin, brg, vcv, covariances):
        self.name = name
        self.description = description
        self.con = con
        self.lat = lat
        self.lon = lon
        self.ohgt = ohgt
        self.ehgt = ehgt
        self.sd_e = sd_e
        self.sd_n = sd_n
        self.sd_u = sd_u
        self.hpu = hpu
        self.vpu = vpu
        self.smaj = smaj
        self.smin = smin
        self.brg = brg
        self.vcv = vcv
        self.covariances = covariances

    def xyz(self):
        return gc.llh2xyz(self.lat, self.lon, self.ehgt)

    def grid(self):
        return gc.geo2grid(self.lat, self.lon)


class DynaMetadata(object):
    def __init__(self, epoch, reference_frame, geoid_model, version):
        self.epoch = epoch
        self.reference_frame = reference_frame
        self.geoid_model = geoid_model
        self.version = version
        # self.file_name = file_name
        # self.file_date = file_date
