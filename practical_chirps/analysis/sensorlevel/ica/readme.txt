Mit preprocessing_batch4all können die MEG-Daten in fif-files vorverarbeitet werden.
Die Hauptskripte sind dabei:
(1) preprocessing_main.m
(2) artifactdetection01.m
(3) artifactdetection02.m
(4) artifactdetection03.m
(5) artifactdetection04.m

(6) helper_function\check_identified_independent_components

Allen Skripten ist gemeinsam, dass oben die Parameter eingestellt werden können. Diese sind von Skript zu Skript verschieden 
und werden im folgenden etwas näher erläutert.

Zu (1):
- Subjektname verweist auf das angelegte m.file mit den Probandeninformation.
- files2preproc ist der Name für die fif-files, die vorverarbeitet werden sollen. Dabei können sowohl einzelne files vorverarbeitet werden
  z.B. Story1Part1 etc. aber auch alle files aufeinmal, so wie es gewünscht ist. Dazu wird Stories angegeben. Für beide gibt es auch die
  Option gemaxfilterte Daten zu laden, dann lautet der Name entsprechende Story1Part1MaxFilter, StoriesMaxFilter. 
- epoch ermöglicht das epochieren der daten in z.B. 10s trials. Tut man das nicht, so liegen die Daten kontinuierlich vor. Verwendet man dabei
  mehrer files auf einmal, so wird jedem file ein trial zugeordnet.
- save_data ermöglicht das Speichern wichtiger erstellter Daten. Die Daten werden dabei immer überschrieben. Diese Option sollte grundsätzlich 
  mit 1 aktiviert sein.
- reject ermöglicht das Entfernen vorher ermittelter ICA-Komponenten.
- sensortype2clean gibt die Sensortypen an bei denen eine ICA-durchgeführt wurde. Bei diesen Sensortypen werden entsprechend die Komponenten gelöscht.
  Zur Verfügung stehen dabei (eeg,megmag,megplanar,meg). Das Einstellen dieser Option sollte sich mit der in artifactdetection01.m durchgeführten ICA
  decken!!!
- Außerdem gibt es die Option während der Durchführung des Skripts nicht benutzte Daten mit c=1 zu entfernen, was möglichweise die Performance erhöht und den
  verwendeten Arbeitsspeicher reduziert.
Zu (2):
- siehe (1)
- ica_sensors gibt die Sensortypen an bei denen eine ICA durchgeführt wird. Dieses sollte sich in (1) mit sensortype2clean decken!!!
- mit check_rank kann der Rank der Datenmatrix berechnet werden. Dieser ist bei gemaxfilterten Daten reduiert (für meg 65 statt 306). 
  Dies führt dazu, dass in diesem Fall die Daten vor der ICA mit einer PCA auf 65 Komponenten reduziert werden.
Zu (3):
- siehe (1)
Zu (4)
- siehe (1)
zu (5)
- siehe (1)
- Mit diesem Skript werden die artefaktkorrelierten ICA-Komponenten gefunden und müssen (noch) manuell notiert werden. Diese werden dann in die Subjectfile manuell
  in die entsprechenen Stellen eingetragen. Da vorgesehen ist, die ICA auf allen Files gleichzeitig durchzuführen, gibt es nur zwei Einträge für die Stories, 
  nämlich einen für gemaxfilterte Daten und einen für nicht-gemaxfilterte Daten. Weiterhin gibt es auch zwei Einträge für das OLSA-file.
- Zum Sichten der Komponenten gibt es folgende Optionen:
  - Kohärenzanalyse zwischen Artefaktkanälen und ICA-Komponenten der Sensortypen (FT der Kreuzkorrelation)
  - Korrelationskoeffizient zwischen Artefaktkanälen und ICA-Komponenten
  - time-locked Analysis (Mittelung über trials). Nur möglich für EKG-Kanal, da alle trials gleichlang sind.
  - Zusätzlich können mit dem Skript check_identified_independent_components.m in helper_functions die Zeitreihen zwischen den gefundenen Komponenten und den 
    Artefaktkanälen angesehen werden. Ebenso kann mit einem Topoplot der Komponenten beobachtet werden, ob die Komponenten typische Muster wie wie Augenartefakten
    und Herzartefakten aufweisen. Damit das funktioniert müssen die gefundenen Komponenten für jeden Sensortyp oben im header des Skripts angegeben werden.
  - Die dann bestätigten Komponenten werden dann in das Subjectfile eingetragen.

zu (6)
- Mit diesem Skript findet die visuelle Inspektion der gefundenen Komponenten statt. Dabei werden für die vorher
  aufgefundenen Komponenten Topoplots und Zeitreihen abgebildet.

Mit den weiteren Skript check_preprocessing.m in helper_functions kann die Wirksamkeit der ICA untersucht werden, indem Daten mit/ohne ICA verglichen werden.
Außerdem kann zusätzlich mit den unverarbeiteten Rohdaten verglichen werden. Zudem können auch die Spektren betrachtet werden.

