%
% make pairs of transient stimuli with gap in between
%

fsamp = 44100;       % Hertz
gap_length = 0.016;   % sec
ngap = round(gap_length*fsamp);
gap = zeros(1,ngap);

s1 = audioread('click.wav');
s2 = audioread('gup30.wav');

% find last non-zero sample in s1 and s2

l1 = length(s1);
for i = l1:-1:1
	if s1(i) ~= 0
		sample1 = i;
		break;
        end;
end;
l2 = length(s2);
for i = l2:-1:1
	if s2(i) ~= 0
		sample2 = i;
		break;
        end;
end;

stimulus = [s1(1:sample1)' gap s2(1:sample2)' gap];
plot(stimulus);

audiowrite('sound1.wav',stimulus, fsamp);

