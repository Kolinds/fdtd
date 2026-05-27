for qTime in range (0, cf.TOTAL_TIME):
    grid.update_Hyfield()
    grid.apply_hyTFSF(incf.ricker, cf.TFSF_BOUNDARY, qTime, 50, 0, 0, cf.STEPS_WAVELENGTH, cf.RICKER_DELAY)
    
    grid.update_Ezfield()
    grid.apply_ezTFSF(incf.ricker, cf.TFSF_BOUNDARY, qTime, 50, 0.5, -0.5, cf.STEPS_WAVELENGTH, cf.RICKER_DELAY)


    grid.abc.second_order()


    grid.r_DFT(qTime)

    hdf5_handler.update_file(qTime, grid.ez, grid.hy)