import requests
import os


def load(atom, charges):
    items = ['ai', 'ce', 'ci', 'en', 'rr', 'tr']

    if not os.path.exists(atom):
        os.mkdir(atom)

    for it in items:
        for ch in charges:
            filename = '{}/{}{}a.{}.tar.gz'.format(atom, atom, ch, it)

            # https://www-amdis.iaea.org/FAC/Si/Si01a.ce.tar.gz
            url = 'https://www-amdis.iaea.org/FAC/{}'.format(filename)
            print(url)
            res = requests.get(url)
            if res.status_code == requests.codes.ok:
                with open(filename, 'wb') as f:
                    f.write(res.content)


if __name__ == '__main__':
    load('He', ['01', '02'])
    load('Li', ['01', '02', '03'])
    load('Be', ['01', '02', '03', '04'])
    load('B', ['01', '02', '03', '04', '05'])
    load('F', ['01', '02', '03', '04', '05', '06'])
    load('C', ['01', '02', '03', '04', '05', '06', '07'])
    load('Ni', ['01', '02', '03', '04', '05', '06', '07', '08'])
    load('O', ['01', '02', '03', '04', '05', '06', '07', '08', '09'])
    load('F', ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10'])
    load('Ne', ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                '11'])
    load('Na', ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                '11', '12'])
    load('Mg', ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                '11', '12', '13'])
    load('Al', ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                '11', '12', '13', '14'])
    load('Si', ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                '11', '12', '13', '14', '15'])
