%% Terminate only at the prescribed observation time.
function stop=stopcheck(timerun,timeend)
stop=timerun(end)>=timeend;
end
