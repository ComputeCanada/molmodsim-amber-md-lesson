import sys
import pandas as pd
import sqlite3
from sqlite3 import Error

print("Usage:")
print("python benchmarks_DB.py \"'AMBER', '20', 'amber/20.9-20.15', 'cedar', 'Broadwell', 'no', 8, 1, 0, 300000, 2.1, NULL, NULL\"\n")

conn = sqlite3.connect('MD_benchmarks.sqlite')
c = conn.cursor()

c.execute("CREATE TABLE IF NOT EXISTS benchmarks (\
id INTEGER PRIMARY KEY AUTOINCREMENT,\
software TEXT NOT NULL,\
version TEXT NOT NULL,\
module TEXT,\
system TEXT NOT NULL,\
cpu_arch TEXT NOT NULL,\
nvlink TEXT NOT NULL,\
ncores INTEGER NOT NULL,\
ntasks INTEGER NOT NULL,\
ngpus INTEGER NOT NULL,\
natoms INTEGER NOT NULL,\
rate REAL NOT NULL,\
command TEXT)")

c.execute("INSERT INTO benchmarks (\
    software,\
    version,\
    module,\
    system,\
    cpu_arch,\
    nvlink,\
    ncores,\
    ntasks,\
    ngpus,\
    natoms,\
    rate,\
    command\
    )VALUES" + "(" + sys.argv[1] + ")") 
conn.commit()


df = pd.read_sql_query("SELECT DISTINCT * from benchmarks ORDER BY rowid DESC LIMIT 1;", conn)
print(df.head())

c.close()
conn.close()
