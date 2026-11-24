def write_global(frame, path):
    import struct
    from compact_io import array, read_n
    with path.open('rb') as src:
        step,np_,nt,t,held,*metrics=struct.unpack('<3i6d',src.read(60))
        if step!=frame.attrs['step'] or t!=frame.attrs['time_s'] or held!=frame['mi'].attrs['held_s']:
            raise ValueError('Global MI clock mismatch')
        group=frame.create_group('global_mi')
        group.attrs.update(held_s=held,longitude='Earth-fixed east longitude, 0..360 exclusive',
                           latitude='geographic latitude, 90..-90 inclusive',phi_gauge='area-weighted full-sphere mean zero',
                           final_residual=metrics[0],removed_mean_A_m2=metrics[1],
                           net_fac_before_A=metrics[2],net_fac_after_A=metrics[3])
        for name in ('potential_V','pedersen_S','hall_S','outward_fac_A_m2'):
            array(group,name,read_n(src,'<f8',(nt,np_)))
        if src.read(1):raise ValueError('Unexpected global MI trailing data')
