%% Stop after the specified observation time.
function stop=stopcheck(timerun,timeend)
stop=timerun(end)>=timeend;
end
