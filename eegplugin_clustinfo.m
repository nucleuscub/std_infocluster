function eegplugin_clustinfo(fig, try_strings, catch_strings)

vers = 'Infocluster_1.0_beta';

% add menu to clustedit
% ---------------------
menu = findobj(fig, 'Label', 'Edit/plot clusters');
structure.uilist = { { } ...
    {'style' 'pushbutton' 'string' 'InfoCluster (beta)' 'Callback' 'args = { STUDY ALLEEG}; STUDY = pop_std_clustinfo(args{:}); clear args;' } { } { } };
structure.geometry = { [1] [1 0.3 1] };
arg = vararg2str( { structure } );
cb_clustedit   = [ trystrs.no_check 'ALLEEGTMP = ALLEEG; [STUDYTMP LASTCOM] = pop_clustedit(STUDY, ALLEEG, [], ' arg ');' catchstrs.update_study];
set(menu, 'callback', cb_clustedit);