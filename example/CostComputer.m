function [cost, cost_state, cost_input] = CostComputer(state,input,options)

cost_state = 0;
cost_input = 0;

for count=1:size(state,2)-1
    cost_state = cost_state+ (state(:,count)-options.xref)'*options.Q*(state(:,count)-options.xref);
end
cost_state = cost_state + (state(:,end) - options.xref)'*options.P*(state(:,end) - options.xref);

for count=1:size(input,2)
    cost_input = cost_input+ input(:,count)'*options.R*input(:,count);
end
  
cost = cost_state + cost_input;

end