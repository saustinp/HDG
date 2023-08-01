Er_loaded = mesh.dgnodes(:,3,:);
Ez_loaded = mesh.dgnodes(:,4,:);

Er_solution = UDG_history(:,8,:,4);
Ez_solution = UDG_history(:,12,:,4);

max(max(Ez_loaded))
max(max(Ez_solution))
max(max(Ez_loaded-Ez_solution))
