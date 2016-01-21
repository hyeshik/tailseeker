import tailseqext
import numpy as np

orig = np.array([[1,2,3,4], [901, 17281, -999, 1]], np.int16)
r = tailseqext.encode_intensity(orig)
print orig
print r
print tailseqext.decode_intensity(r)



decoder = tailseqext.IndexDecoder(['ACCGGT', 'GATTCA', 'CAGATC'])
print decoder('ACCGAT')
print decoder('GATTCA')
print decoder('GATTCC')
#print decoder('ACCTCA')



loc = tailseqext.PolyALocator({'A': -10, 'G': -10, 'C': -10, 'N': -5, 'T': 1})
print loc('GTTTTTTTTTTTAACGTAA', 4)
