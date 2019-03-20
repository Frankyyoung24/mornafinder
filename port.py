#!/usr/bin/env python

from __future__ import print_function

import math
import os
import re
import sys
import time
from string import maketrans

from colorama import Fore, Style, init

# Colorama
init()


def localtime(ctime=None):
    t = time.localtime(ctime)
    return (t.tm_sec, t.tm_min, t.tm_hour, t.tm_mday, t.tm_mon - 1,
            t.tm_year - 1900, t.tm_wday + 1, t.tm_yday - 1, t.tm_isdst)


def substr(string, beg, length=None):
    '''
    Get a substring from string
    @param string    The string to be cut
    @param beg       The start offset
    @param length    The substring length from beg offset to be cut
    '''
    assert isinstance(string, str) or isinstance(
        string, unicode), 'string must be str/unicode type'
    if length is not None:
        if length < 0:
            return string[beg:length]
        else:
            return string[beg:beg + length]
    else:
        return string[beg:]


def pprint(*string):
    '''
    Print without newline at the end
    '''
    print(*string, end='')


def hash_sort_key(hsh, func, **kwargs):
    '''
    Sort a hash by func, and get the sorted keys
    '''
    assert isinstance(hsh, dict), 'hsh must be a dict'
    return [x[0] for x in sorted(hsh.items(), key=func, **kwargs)]


def hash_defined(hsh, *keys):
    '''
    Check a dict key chain value exist or not.
    (Only None value is excluded, '', 0, False will be treated as defined)
    eg,
    hsh = {}, hash_value_true(hsh, 'a', 'b') => False
    hsh = {'a': {'b': 1}}, hash_value_true(hsh, 'a', 'b') => True
    '''
    if isinstance(hsh, dict):
        tmp = hsh
        for key in keys:

            if not (isinstance(tmp, dict) and key in tmp.keys() and tmp[key] is not None):
                return False
            else:
                tmp = tmp[key]
        return True
    else:
        False


def hash_value_true(hsh, *keys):
    '''
    Check a dict key chain value exist, and have a True value or not.
    eg,
    hsh = {}, hash_value_true(hsh, 'a', 'b') => False
    hsh = {'a': {'b': 1}}, hash_value_true(hsh, 'a', 'b') => True
    '''
    if isinstance(hsh, dict):
        tmp = hsh
        for key in keys:
            if not (isinstance(tmp, dict) and key in tmp.keys() and tmp[key]):
                return False
            else:
                tmp = tmp[key]
        return True
    else:
        False


def create_hash_key_chain(hsh, endvalue, *keys):
    '''
    Create a perl-like hash chain with specified keys, if hash value not exist then create one.
    @param hsh         the hash(dict) value
    @param endvalue    the value that to be assigned to the last node
    @param *keys       a list for keys
    eg, hash = {}, you need to use this stmt: ${hash}{'a'}{'b'}{'c'} += 1,
    then you should call create_hash_key_chain(hash, 0, 'a', 'b', 'c') to init the node value
    '''
    assert isinstance(hsh, dict), 'hsh must be dict type'
    temp = hsh
    added = False
    for i in range(len(keys)):
        if isinstance(temp, dict):
            if keys[i] in temp.keys():
                added = False
            else:
                temp[keys[i]] = {}
                added = True

        if i == len(keys) - 1:  # last key
            if added:
                temp[keys[i]] = endvalue
        else:
            temp = temp[keys[i]]


def fileparse(filename, pattern=None):
    result = os.path.split(filename)
    if not pattern:
        return (result[1], result[0], '')
    else:
        name = re.split(pattern, result[1])[0]
        m = re.findall(pattern, result[1])
        ext = ''
        if m:
            ext = m[0]
        return (name, result[0], ext)


def chop(string):
    '''
    Remove the last character from string
    '''
    assert isinstance(string, str) or isinstance(
        string, unicode), 'string must be str/unicode type'
    string = string[:-1]
    return string


def str_reverse(string):
    '''
    Reverse the given string. eg, str_reverse('123') => '321'
    '''
    arr = list(string)
    arr.reverse()
    return ''.join(arr)


def file_s(filePath):
    '''
    if -s xx
    '''
    assert isinstance(filePath, str) or isinstance(
        filePath, unicode), 'filePath should be str/unicode'
    if os.path.exists(filePath) and os.path.getsize(filePath) > 0:
        return True
    else:
        return False


def esplit(data):
    '''
    split()
    '''
    arr = re.split(r'\s', data)
    return filter(lambda x: x != '', arr)


def ssplit(data):
    '''
    split(//)
    '''
    return list(data)


def print_stderr(*messages):
    '''
    Print messags to stderr
    '''
    print(*messages, sep='', file=sys.stderr, end='')


def die(*messages):
    '''
    Print a message, and exit the script
    '''
    print(*messages, sep='', end='')
    sys.exit(-1)


def open_or_die2(filename, mode):
    '''
    Open the file with filename, in mode, if IOError occured, then exit
    '''
    if re.search('w', mode):
        return open_or_die(filename, mode, 'can not create file {}\n'.format(filename))
    else:
        return open_or_die(filename, mode, 'can not open file {}\n'.format(filename))


def open_or_die(filename, mode, die_message):
    '''
    Open the file with filename, in mode, if IOError occured, then
    exit with die_message
    '''
    try:
        fh = open(filename, mode)
    except IOError:
        die(die_message)

    return fh


def print_stderr_red(message):
    print_stderr(Fore.RED + message)
    print_stderr(Style.RESET_ALL)


def tr(string, pattern, replacement):
    assert isinstance(string, str), 'string should be str type'
    assert isinstance(pattern, str), 'pattern should be str type'
    assert isinstance(replacement, str), 'replacement should be str type'
    # assert len(pattern) <= len(replacement), 'pattern length should less than replacement'

    if replacement:
        if len(replacement) < len(pattern):
            temp = replacement[-1] * (len(pattern) - len(replacement))
            replacement += temp
        elif len(replacement) > len(pattern):
            replacement = replacement[:len(pattern)]

        tran_tbl = maketrans(pattern, replacement[:len(pattern)])
        return string.translate(tran_tbl)
    else:
        temp = string
        for i in list(pattern):
            temp = temp.replace(i, '')
        return temp


def Nicenumber(num):
    numarr = re.split(r'[.,]', str(num))
    number = numarr[0]
    n = ssplit(str_reverse(number))
    res = ''
    for j in range(0, len(number)):
        if j % 3 == 0 and j != 0:
            res += '.{}'.format(n[j])
        else:
            res += '{}'.format(n[j])

    res = str_reverse(res)
    if len(numarr) > 1 and numarr[1]:
        res += ',{}'.format(numarr[1])

    return res


def mean(arr):
    '''
    Calculate the mean of values in arr
    @param    arr - list/number
    @retval   The mean
    '''
    if not(isinstance(arr, tuple) or isinstance(arr, list) or
           isinstance(arr, int) or isinstance(arr, float) or
           isinstance(arr, long)) or not arr:
        return None

    if isinstance(arr, int) or isinstance(arr, long) or isinstance(arr, float):
        return arr

    sum = 0

    for item in arr:
        sum += item

    return (float(sum) / len(arr))


def sd(arr):
    '''
    Calculate the standard deviation from given number of array
    @param    arr - list/tuple/number
    @retval   Standard deviation value
    '''
    if not(isinstance(arr, tuple) or isinstance(arr, list) or
           isinstance(arr, int) or isinstance(arr, float) or
           isinstance(arr, long)) or not arr:
        return None

    if isinstance(arr, int) or isinstance(arr, long) or isinstance(arr, float):
        return arr

    m = mean(arr)
    sum_squares = 0
    for var in arr:
        dif = m - var
        square = dif ** 2
        sum_squares += square

    freedom = len(arr) - 1
    _sd = 0
    if freedom != 0:
        _sd = math.sqrt(float(sum_squares) / freedom)

    return _sd


def pround(num):
    '''
    Round a float to integer
    '''
    return int(num + 0.5)


def round_decimals(num, dec):
    '''
    Round a decimal number to given decimal position
    '''
    if dec == 0:
        return pround(num)
    else:
        return pround(10**dec * num) / float(10**dec)


def div(num, den):
    '''
    Number division
    '''
    if den == 0:
        return 0
    else:
        return float(num) / den


def max2(a, b):
    '''
    Get the maximum one upon given two numbers
    '''
    return a if a > b else b


def min2(a, b):
    '''
    Get the minimum one upon given two numbers
    '''
    return b if a > b else a


def rev(seq):
    '''
    reverses the order of nucleotides in a sequence
    '''
    return str_reverse(seq)


def com(seq):
    '''
    the complementary of a sequence
    '''
    return tr(seq, 'acgtuACGTU', 'TGCAATGCAA')


def revcom(seq):
    '''
    reverse complement
    '''
    return rev(com(seq))


if __name__ == '__main__':
    import unittest

    class PortTest(unittest.TestCase):

        @classmethod
        def tearDownClass(cls):
            os.system('rm test.txt')

        def test_str_reverse(self):
            string = 'abcdefg'
            stringR = str_reverse(string)
            self.assertEqual(stringR, 'gfedcba')

        def test_ssplit(self):
            m = 'who am i\nhello world'
            self.assertListEqual(ssplit(m), ['w', 'h', 'o', ' ', 'a', 'm', ' ',
                                             'i', '\n', 'h', 'e', 'l', 'l', 'o',
                                             ' ', 'w', 'o', 'r', 'l', 'd'])

        def test_esplit(self):
            m = ''' who\tam i\r\n i am enix'''
            arr = esplit(m)
            self.assertEqual(len(arr), 6)
            self.assertListEqual(arr, ['who', 'am', 'i', 'i', 'am', 'enix'])

            m = 'who'
            arr = esplit(m)
            self.assertEqual(len(arr), 1)
            self.assertListEqual(arr, ['who'])

        def test_chop(self):
            m = '123'
            self.assertEqual(chop(m), '12')
            m = ''
            self.assertEqual(chop(m), '')

        def test_hash_sort_key(self):
            d = {'a': {'k': 8}, 'b': {'k': 4}, 'c': {'k': 0}}
            self.assertListEqual(hash_sort_key(
                d, lambda x: x[1]['k']), ['c', 'b', 'a'])
            d = {'a': {'k': 8}, 'b': {'k': 4}, 'c': {'k': 0}}
            self.assertListEqual(hash_sort_key(
                d, lambda x: x[1]['k'], reverse=True), ['a', 'b', 'c'])
            d = {'a': {'k': 8}}
            self.assertListEqual(hash_sort_key(d, lambda x: x[1]['k']), ['a'])
            d = {}
            self.assertListEqual(hash_sort_key(d, lambda x: x[1]['k']), [])

        def test_create_hash_key_chain(self):
            hsh = {}
            create_hash_key_chain(hsh, 0, 'a', 'b', 'c')
            self.assertDictEqual(hsh, {'a': {'b': {'c': 0}}})

            hsh = {}
            create_hash_key_chain(hsh, 0, 'a')
            self.assertDictEqual(hsh, {'a': 0})

            # if child node in the middle of chain is not exist
            hsh = {'a': {'b': {}}}
            create_hash_key_chain(hsh, 0, 'a', 'b', 'c')
            self.assertDictEqual(hsh, {'a': {'b': {'c': 0}}})

            hsh = {'a': {'b': {}}}
            create_hash_key_chain(hsh, None, 'a', 'b', 'c')
            self.assertDictEqual(hsh, {'a': {'b': {'c': None}}})

            # if child node in the middle of chain is not dict, then no changes
            # is applied
            hsh = {'a': {'b': 1}}
            create_hash_key_chain(hsh, 1000, 'a', 'b', 'c')
            self.assertDictEqual(hsh, {'a': {'b': 1}})

            # if node exist, then endvalue will not assigned
            hsh = {'a': {'b': {'c': 1}}}
            create_hash_key_chain(hsh, 999, 'a', 'b', 'c')
            self.assertDictEqual(hsh, {'a': {'b': {'c': 1}}})

        def test_nice_number(self):
            num = '1234'
            self.assertEqual(Nicenumber(num), '1.234')
            num = '123'
            self.assertEqual(Nicenumber(num), '123')
            num = '123,41'
            self.assertEqual(Nicenumber(num), '123,41')
            num = '1234,56'
            self.assertEqual(Nicenumber(num), '1.234,56')
            num = '1234567890'
            self.assertEqual(Nicenumber(num), '1.234.567.890')

        def test_mean(self):
            a = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            self.assertEqual(mean(a), 5.5)
            a = []
            self.assertEqual(mean(a), None)
            a = 10
            self.assertEqual(mean(a), a)
            a = None
            self.assertEqual(mean(a), None)
            a = ''
            self.assertEqual(mean(a), None)
            a = {}
            self.assertEqual(mean(a), None)

        def test_sd(self):
            a = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            self.assertEqual(sd(a), 3.0276503540974917)
            a = 1
            self.assertEqual(sd(a), a)
            a = None
            self.assertEqual(sd(a), None)
            a = False
            self.assertEqual(sd(a), None)
            a = ''
            self.assertEqual(sd(a), None)
            a = []
            self.assertEqual(sd(a), None)
            a = {}
            self.assertEqual(sd(a), None)

        def test_div(self):
            self.assertEqual(div(10, 0), 0)
            self.assertEqual(div(10, 5), 2)
            self.assertEqual(div(10, 3), float(10) / 3)

        def test_max2(self):
            self.assertEqual(max2(10, 5), 10)
            self.assertEqual(max2(10.12313, 10.321), 10.321)
            self.assertEqual(max2(10, 0), 10)

        def test_min2(self):
            self.assertEqual(min2(10, 5), 5)
            self.assertEqual(min2(10.123, 10.321), 10.123)
            self.assertEqual(min2(10, 0), 0)

        def test_pround(self):
            a = 1
            self.assertEqual(pround(a), a)
            a = 10.3
            self.assertEqual(pround(a), 10)
            a = 10.5
            self.assertEqual(pround(a), 11)

        def test_round_decimals(self):
            a = 10.1458910
            self.assertEqual(round_decimals(a, 2), 10.15)
            a = 10.512761
            self.assertEqual(round_decimals(a, 2), 10.51)
            self.assertEqual(round_decimals(a, 1), 10.5)
            self.assertEqual(round_decimals(a, 0), 11)

        def test_tr(self):
            m = 'xyz12345980.\\xyz'
            mm = tr(m, '[x1.]', '[1x0]')
            self.assertEqual(mm, '1yzx23459800\\1yz')

            mm = tr(m, '[x1.]', '[1x0]321321321')
            self.assertEqual(mm, '1yzx23459800\\1yz')

            mm = tr(m, '[x1.\\]', '[1x0.]')
            self.assertEqual(mm, '1yzx23459800.1yz')

            mm = tr(m, 'x', '')
            self.assertEqual(mm, 'yz12345980.\\yz')

            mm = tr(m, 'xyz', '')
            self.assertEqual(mm, '12345980.\\')

            mm = tr(m, 'xyz', '1')
            self.assertEqual(mm, '11112345980.\\111')

        def test_open_or_die_should_ok(self):
            m = 'cannot create test.txt'
            fh = open_or_die('test.txt', 'w+', m)
            fh.write('hello')
            fh.close()
            self.assertTrue(True)

        def test_hash_defined(self):
            h = {}
            self.assertFalse(hash_defined(h, 'a', 'b'))

            h = {'a': None}
            self.assertFalse(hash_defined(h, 'a'))

            h = {'a': 0}
            self.assertTrue(hash_defined(h, 'a'))

            h = {'a': False}
            self.assertTrue(hash_defined(h, 'a'))

            h = {'a': ''}
            self.assertTrue(hash_defined(h, 'a'))

            h = {'a': {}}
            self.assertTrue(hash_defined(h, 'a'))

            h = {'a': []}
            self.assertTrue(hash_defined(h, 'a'))

            h = {'a': {'b': 1}}
            self.assertTrue(hash_defined(h, 'a', 'b'))

            h = {'a': {'b': None}}
            self.assertFalse(hash_defined(h, 'a', 'b'))

            h = None
            self.assertFalse(hash_defined(h))

        def test_substr(self):
            a = 'abcdefghijk'
            self.assertEqual(substr(a, -8, 4), 'defg')
            self.assertEqual(substr(a, -8, 40), 'defghijk')
            self.assertEqual(substr(a, 3, 4), 'defg')
            self.assertEqual(substr(a, 0, 4), 'abcd')
            self.assertEqual(substr(a, 3, 0), '')
            self.assertEqual(substr(a, -5, -2), 'ghi')
            self.assertEqual(substr(a, -5), 'ghijk')

        def test_hash_value_true(self):
            h = {}
            self.assertFalse(hash_value_true(h, 'a', 'b'))

            h = {'a': 1}
            self.assertTrue(hash_value_true(h, 'a'))

            h = {'a': {}}
            self.assertFalse(hash_value_true(h, 'a', 'b'))

            h = {'a': {'b': 1}}
            self.assertTrue(hash_value_true(h, 'a', 'b'))

            h = {'a': {'b': None}}
            self.assertFalse(hash_value_true(h, 'a', 'b'))

            h = {'a': {'b': {}}}
            self.assertFalse(hash_value_true(h, 'a', 'b'))

            h = 1
            self.assertFalse(hash_value_true(h, 'a', 'b'))

            h = None
            self.assertFalse(hash_value_true(h, 'a', 'b'))

            h = 'abc'
            self.assertFalse(hash_value_true(h, 'a', 'b'))

        def test_localtime(self):
            (second, minute, hour, dayOfMonth, month, yearOffset,
             dayOfWeek, dayOfYear, daylightSavings) = localtime()
            # self.assertEqual(hour, 14)
            # self.assertEqual(dayOfMonth, 28)
            # self.assertEqual(month, 2)
            # self.assertEqual(yearOffset, 116)
            # self.assertEqual(dayOfWeek, 1),
            # self.assertEqual(dayOfYear, 87)
            # self.assertEqual(daylightSavings, 0)

        def test_pprint(self):
            pprint('Test', 'pprint')
            pprint('\n')

    unittest.main()
