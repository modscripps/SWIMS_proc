function DisplayProgress(counter,freq)
%function DisplayProgress(counter,freq)
%Every (freq) times display a note saying that that number is being processed.
%Pass the counter to be checked and the frequency with which to do so.
if rem(counter,freq)==0
	disp (['Working on number ' num2str(counter)])
end
