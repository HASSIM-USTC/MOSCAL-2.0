import subprocess

num_of_threads = 7

magic_word_list = [
    'bose-dt', 'bose-corr-dt', 'bose-corr-dt-f',
    'bsm-corr-dt-filter',
    'bose-quad-cfw-dt-filter', 'bose-quad-cfw-backward-dt-filter', 'bose-quad-corr-dt-filter',
    'bose-quad-corr-dt-filter-f', 'bose-quad-dt-filter', 'bose-quad-imag-dt-filter',
    'fermi-corr-dt-filter', 'fermi-corr-dt-hei-filter',
    'fermi-corr-dt-hei-noise-filter', 'fermi-corr-dt-noise-filter',
    'fermi-corr-dw-filter', 'fermi-corr-dw-hei-filter',
    'fermi-corr-dw-filter-no-sym', 'fermi-corr-dw-hei-filter-no-sym',
    'fermi-corr-dw-hei-noise-filter', 'fermi-corr-dw-noise-filter',
    'fermi-dt-filter', 'fermi-sc2-filter', 'fermi-dt-imag-filter',
    'fermi-dt-cfw-filter', 'fermi-corr-dt', 'fermi-corr-dw', 'fermi-corr-dt-noise', 'fermi-corr-dw-noise', 'fermi-sc2', 'fermi-dt',
    'fermi-quad-dt', 'fermi-quad-corr-dt',
]

print(len(magic_word_list))

for magic_word in magic_word_list:
    print(magic_word)
    cmd = "python3 -W ignore {}.py {}".format(
        magic_word, num_of_threads)
    result = subprocess.call(cmd, shell=True)
    print("\n")

cmd = "rm time-*"
cmd += '\n' + "rm prop-*"
cmd += '\n' + "rm curr.dat"
cmd += '\n' + "rm out-*"
cmd += '\n' + "rm input.json"
subprocess.call(cmd, shell=True)
