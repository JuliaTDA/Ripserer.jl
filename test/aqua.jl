using Aqua
using TDARipserer
using TDAPersistenceDiagrams

Aqua.test_all(TDARipserer; ambiguities=false, piracies=false)
Aqua.test_piracies(TDARipserer; treat_as_own=[PersistenceDiagram, PersistenceInterval])
