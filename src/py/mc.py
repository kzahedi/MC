# This file is part of Quantifying Morphological Computation (MC)
# License is MIT, see LICENCE and http://github.com/kzahedi/MC/LICENSE

###########################################################################
#                            binning functions                            #
###########################################################################

def bin_value(_v, _bins):
  return min(int(_v * _bins), _bins-1)

# TODO
def bin_vector(_v, _bins):
  return [min(int(_v[i] * _bins), _bins-1) for i in range(0,3)]

# TODO
def bin_scaled_data_for_each_marker(_data, _bins):
  _new_data = {}
  for _key in _data.keys():
    _new_data[_key] = [bin_vector(_v, _bins) for _v in _data[_key]]
  return _new_data

def combine_bin_vector(_v, _bins):
    return sum([_v[i] * pow(_bins,i) for i in range(0, len(_v))])

def unique_valued_list(_lst):
    _myset = set(_lst)
    return list(_myset)

def relabel_vector(_lst):
    _mylst = unique_valued_list(_lst)
    return [_mylst.index(_v) for _v in _lst]

def combine_bins_for_each_marker(_data, _bins):
    _new_data = {}
    for _key in _data.keys():
        _new_data[_key] = relabel_vector([combine_bin_vector(_v, _bins) for _v in _data[_key]])
    return _new_data

def combine_random_variables(_lst_of_lsts, _bins):
    return relabel_vector([combine_bin_vector([_v[i] for _v in _lst_of_lsts], _bins) for i in range(0, len(_lst_of_lsts[1]))])

###########################################################################
#                   calculating probabilities from data                   #
###########################################################################

def emperical_joint_distribution(_w_prime, _w, _a):
  _p = numpy.zeros((max(_w_prime)+1, max(_w)+1, max(_a)+1))

  _l = len(_w_prime)
  for index in range(0, _l):
    _p[_w_prime[index], _w[index], _a[index]] = _p[_w_prime[index], _w[index], _a[index]] + 1.0

  for _i_index in range(0, _p.shape[0]):
    for _j_index in range(0, _p.shape[1]):
      for _k_index in range(0, _p.shape[2]):
        _p[_i_index,_j_index,_k_index] = _p[_i_index,_j_index,_k_index] / float(_l)

  _s = sum(sum(sum(_p)))
  _p = _p / _s
  return _p

def calc_p_w_prime_given_w(_joint_distribution):
    _p_w_prime_w = _joint_distribution.sum(axis=2)
    _p_w         = _joint_distribution.sum(axis=(0,2))
    for _w_prime in range(0,_joint_distribution.shape[0]):
        for _w in range(0, _joint_distribution.shape[1]):
            _p_w_prime_w[_w_prime, _w] = _p_w_prime_w[_w_prime, _w] / _p_w[_w]
    return _p_w_prime_w

def calc_p_w_prime_given_a(_joint_distribution):
    _p_w_prime_a = _joint_distribution.sum(axis=1)
    _p_a         = _joint_distribution.sum(axis=(0,1))
    for _w_prime in range(0,_joint_distribution.shape[0]):
        for _a in range(0, _joint_distribution.shape[2]):
            if _p_w_prime_a[_w_prime, _a] != 0.0 and _p_a[_a] != 0.0:
                _p_w_prime_a[_w_prime, _a] = _p_w_prime_a[_w_prime, _a] / _p_a[_a]
    return _p_w_prime_a

def calc_p_w_prime_given_w_a(_joint_distribution):
    _p_w_a               = _joint_distribution.sum(axis=0)
    _p_w_prime_given_w_a = numpy.zeros(_joint_distribution.shape)
    for _w_prime in range(0, _joint_distribution.shape[0]):
        for _w in range(0, _joint_distribution.shape[1]):
            for _a in range(0, _joint_distribution.shape[2]):
                if _joint_distribution[_w_prime, _w, _a] != 0.0 and _p_w_a[_w,_a] != 0.0:
                    _p_w_prime_given_w_a[_w_prime, _w, _a] = _joint_distribution[_w_prime, _w, _a] / _p_w_a[_w,_a]
    return _p_w_prime_given_w_a

###########################################################################
#                           MC quantifications                            #
###########################################################################

def calculate_concept_one(_joint_distribution):
    _p_w_prime_given_w   = calc_p_w_prime_given_w(_joint_distribution)
    _p_w_prime_given_w_a = calc_p_w_prime_given_w_a(_joint_distribution)
    _r = 0
    for _w_prime in range(0, _joint_distribution.shape[0]):
        for _w in range(0, _joint_distribution.shape[1]):
            for _a in range(0, _joint_distribution.shape[2]):
                if _joint_distribution[_w_prime, _w, _a] != 0.0 and \
                   _p_w_prime_given_w[_w_prime, _w]      != 0.0 and \
                   _p_w_prime_given_w_a[_w_prime, _w, _a] != 0.0:
                    _r = _joint_distribution[_w_prime, _w, _a] * \
                         (math.log(_p_w_prime_given_w_a[_w_prime, _w, _a], 2) -
                          math.log(_p_w_prime_given_w[_w_prime, _w], 2))
    return _r

def calculate_concept_two(_joint_distribution):
    _p_w_prime_given_a   = calc_p_w_prime_given_a(_joint_distribution)
    _p_w_prime_given_w_a = calc_p_w_prime_given_w_a(_joint_distribution)
    _r = 0
    for _w_prime in range(0, _joint_distribution.shape[0]):
        for _w in range(0, _joint_distribution.shape[1]):
            for _a in range(0, _joint_distribution.shape[2]):
                if _joint_distribution[_w_prime, _w, _a] != 0.0 and _p_w_prime_given_a[_w_prime, _a] != 0.0 and _p_w_prime_given_w_a[_w_prime, _w, _a] != 0.0:
                    _r = _joint_distribution[_w_prime, _w, _a] * \
                         (math.log(_p_w_prime_given_w_a[_w_prime, _w, _a], 2) -
                          math.log(_p_w_prime_given_a[_w_prime, _a], 2))
    return _r

# TODO two local variants of MC
