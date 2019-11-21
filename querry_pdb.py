import urllib.request

url = 'http://www.rcsb.org/pdb/rest/search'

queryText = """
<orgPdbQuery>
  <queryType>org.pdb.query.simple.UniprotGeneNameQuery</queryType>
    <query>PRNP</query>
</orgPdbQuery>


"""

print("query:\n", queryText)
print("querying PDB...\n")

req = urllib.request.Request(url, data=queryText.encode("utf-8"))
f = urllib.request.urlopen(req)


result = f.read().decode()
aa = result.split()
aa.sort()

if result:
	for a in aa:
		print(a)
	print("Found number of PDB entries:", result.count('\n'))
else:
	print("Failed to retrieve results" )
    
