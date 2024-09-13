# This file is part of minMutFinder.
#
# minMutFinder is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# minMutFinder is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with minMutFinder. If not, see <https://www.gnu.org/licenses/>.
#
# Copyright (C) 2024 Ignasi Prats Méndez

#!/usr/bin/env python3
import pandas as pd
import os


def grep(file, pattern, tpe):
    aux = False
    result = ''
    if 'inverse' not in tpe:
        if type(file) is str:
            if os.path.isfile(file):
                with open(file, 'r') as infile:
                    for line in infile:
                        if pattern in line:
                            result += line
                            aux = True
                infile.close()
            else:
                for line in file.split('\n'):
                    if pattern in line:
                        result += '\n' + line
                        aux = True

        elif type(file) is dict:
            for i in file:
                if pattern in i:
                    result += i
                    aux = True
        elif type(file) is pd.DataFrame:
            result = pd.DataFrame(columns=file.columns)
            for column in range(0, file.shape[1]-1):
                for row in range(0, file.shape[0]-1):
                    if pattern in file[column][row]:
                        if result.size > 0:
                            result = pd.concat([result, file.iloc[[row]]])
                        else:
                            result = file.iloc[[row]]
                        aux = True
            result = result.reset_index(drop=True)
        else:
            print('Not supproted data type')
        if tpe == 'bool':
            return aux
        elif tpe == 'text':
            return result
        else:
            return result, aux
    else:
        if type(file) is str:
            with open(file, 'r') as infile:
                for line in infile:
                    if pattern not in line:
                        result += line
                        aux = True
            infile.close()
        elif type(file) is dict:
            for i in file:
                if pattern not in i:
                    result += i
                    aux = True
        elif type(file) is pd.DataFrame:
            result = pd.DataFrame(columns=file.columns)
            for column in range(0, file.shape[1]-1):
                for row in range(0, file.shape[0]-1):
                    if pattern not in file[column][row].to_string():
                        if result.size > 0:
                            result = pd.concat([result, file.iloc[[row]]])
                        else:
                            result = file.iloc[[row]]
                        aux = True
            result = result.reset_index(drop=True)
        else:
            print('Not supproted sadfdata type')
        if tpe == 'inverse-bool':
            return aux
        elif tpe == 'inverse-text':
            return result
        else:
            return result, aux
