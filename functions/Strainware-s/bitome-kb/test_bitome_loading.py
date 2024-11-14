from pathlib import Path
from bitome.core import Bitome

# load a bitome object with data and matrix; save it to a pickle; load it again
test_bitome = Bitome(Path('data', 'NC_000913.3.gb'))
test_bitome.load_data(regulon_db=True)
# test_bitome.load_matrix()
# test_bitome.load_pre_compressed_matrix()
test_bitome.save()
test_bitome = Bitome.init_from_file(Path('bitome.pkl'))
# print(test_bitome.validate_matrix())
