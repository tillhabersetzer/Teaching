function questions(story,part)
%--------------------------------------------------------------------------
Box.Interpreter = 'tex';
Box.WindowStyle = 'modal';

%--------------------------------------------------------------------------
if (story==1 && part==1)
 
    msgbox({'\bf\fontsize{36}1.) Was stört den Erzähler am alten Mann?','','','2.) Wie viele Nächte in Folge schaut der Erzähler','     zum alten Mann ins Zimmer?', ...
            '','','3.) Was ist in der achten Nacht anders als in den Nächten davor?'},'Fragenpaket 1',Box);
       
end

%--------------------------------------------------------------------------
if (story==1 && part==2)

    msgbox({'\bf\fontsize{36}1.) Wie wurde der alte Mann umgebracht?','','','2.) Wo befindet sich der Leichnam des alten Mannes?', ...
            '','','3.) Warum gesteht der Erzähler den  Polizisten den Mord?'},'Fragenpaket 2',Box)

end

%--------------------------------------------------------------------------
if (story==2 && part==1)

    msgbox({'\bf\fontsize{36}1.) Was ist mit dem Begriff des Roten Todes gemeint?','','','2.) Was haben der Prinz und sein Gefolge gemacht,','     um dem Roten Tod zu entgehen?', ...
            '','','3.) Wie sieht das siebte Gemach aus und was befindet sich dort?'},'Fragenpaket 1',Box)

end

%--------------------------------------------------------------------------
if (story==2 && part==2)

    msgbox({'\bf\fontsize{36}1.) Wie oft schlägt die Uhr als die mysteriöse Gestalt auf','     dem Maskenball auftaucht?','','','2.) Wie sieht die mysteriöse Gestalt aus?', ...
            '','','3.) Was geschieht zwischen Prinz Prospero und','     der mysteriösen Gestalt?'},'Fragenpaket 2',Box)

end


% oder entwerfe texte als png bild und lade sie mit imshow als foto hoch...

end