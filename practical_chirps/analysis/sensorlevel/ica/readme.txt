Mit preprocessing_batch4all k�nnen die MEG-Daten in fif-files vorverarbeitet werden.
Die Hauptskripte sind dabei:
(1) preprocessing_main.m
(2) artifactdetection01.m
(3) artifactdetection02.m
(4) artifactdetection03.m
(5) artifactdetection04.m

(6) helper_function\check_identified_independent_components

Allen Skripten ist gemeinsam, dass oben die Parameter eingestellt werden k�nnen. Diese sind von Skript zu Skript verschieden 
und werden im folgenden etwas n�her erl�utert.

Zu (1):
- Subjektname verweist auf das angelegte m.file mit den Probandeninformation.
- files2preproc ist der Name f�r die fif-files, die vorverarbeitet werden sollen. Dabei k�nnen sowohl einzelne files vorverarbeitet werden
  z.B. Story1Part1 etc. aber auch alle files aufeinmal, so wie es gew�nscht ist. Dazu wird Stories angegeben. F�r beide gibt es auch die
  Option gemaxfilterte Daten zu laden, dann lautet der Name entsprechende Story1Part1MaxFilter, StoriesMaxFilter. 
- epoch erm�glicht das epochieren der daten in z.B. 10s trials. Tut man das nicht, so liegen die Daten kontinuierlich vor. Verwendet man dabei
  mehrer files auf einmal, so wird jedem file ein trial zugeordnet.
- save_data erm�glicht das Speichern wichtiger erstellter Daten. Die Daten werden dabei immer �berschrieben. Diese Option sollte grunds�tzlich 
  mit 1 aktiviert sein.
- reject erm�glicht das Entfernen vorher ermittelter ICA-Komponenten.
- sensortype2clean gibt die Sensortypen an bei denen eine ICA-durchgef�hrt wurde. Bei diesen Sensortypen werden entsprechend die Komponenten gel�scht.
  Zur Verf�gung stehen dabei (eeg,megmag,megplanar,meg). Das Einstellen dieser Option sollte sich mit der in artifactdetection01.m durchgef�hrten ICA
  decken!!!
- Au�erdem gibt es die Option w�hrend der Durchf�hrung des Skripts nicht benutzte Daten mit c=1 zu entfernen, was m�glichweise die Performance erh�ht und den
  verwendeten Arbeitsspeicher reduziert.
Zu (2):
- siehe (1)
- ica_sensors gibt die Sensortypen an bei denen eine ICA durchgef�hrt wird. Dieses sollte sich in (1) mit sensortype2clean decken!!!
- mit check_rank kann der Rank der Datenmatrix berechnet werden. Dieser ist bei gemaxfilterten Daten reduiert (f�r meg 65 statt 306). 
  Dies f�hrt dazu, dass in diesem Fall die Daten vor der ICA mit einer PCA auf 65 Komponenten reduziert werden.
Zu (3):
- siehe (1)
Zu (4)
- siehe (1)
zu (5)
- siehe (1)
- Mit diesem Skript werden die artefaktkorrelierten ICA-Komponenten gefunden und m�ssen (noch) manuell notiert werden. Diese werden dann in die Subjectfile manuell
  in die entsprechenen Stellen eingetragen. Da vorgesehen ist, die ICA auf allen Files gleichzeitig durchzuf�hren, gibt es nur zwei Eintr�ge f�r die Stories, 
  n�mlich einen f�r gemaxfilterte Daten und einen f�r nicht-gemaxfilterte Daten. Weiterhin gibt es auch zwei Eintr�ge f�r das OLSA-file.
- Zum Sichten der Komponenten gibt es folgende Optionen:
  - Koh�renzanalyse zwischen Artefaktkan�len und ICA-Komponenten der Sensortypen (FT der Kreuzkorrelation)
  - Korrelationskoeffizient zwischen Artefaktkan�len und ICA-Komponenten
  - time-locked Analysis (Mittelung �ber trials). Nur m�glich f�r EKG-Kanal, da alle trials gleichlang sind.
  - Zus�tzlich k�nnen mit dem Skript check_identified_independent_components.m in helper_functions die Zeitreihen zwischen den gefundenen Komponenten und den 
    Artefaktkan�len angesehen werden. Ebenso kann mit einem Topoplot der Komponenten beobachtet werden, ob die Komponenten typische Muster wie wie Augenartefakten
    und Herzartefakten aufweisen. Damit das funktioniert m�ssen die gefundenen Komponenten f�r jeden Sensortyp oben im header des Skripts angegeben werden.
  - Die dann best�tigten Komponenten werden dann in das Subjectfile eingetragen.

zu (6)
- Mit diesem Skript findet die visuelle Inspektion der gefundenen Komponenten statt. Dabei werden f�r die vorher
  aufgefundenen Komponenten Topoplots und Zeitreihen abgebildet.

Mit den weiteren Skript check_preprocessing.m in helper_functions kann die Wirksamkeit der ICA untersucht werden, indem Daten mit/ohne ICA verglichen werden.
Au�erdem kann zus�tzlich mit den unverarbeiteten Rohdaten verglichen werden. Zudem k�nnen auch die Spektren betrachtet werden.

