#!/bin/bash
cat states/$1/procdata/{reactant,saddle,product}_$2.con > temp_mov.con
xyz temp_mov.con
