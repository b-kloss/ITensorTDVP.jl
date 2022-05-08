module params

export outf, tebd_cutoff, tebd_dt,tebd_order, log_dt, final_time
export omega,t, gamma, N, boson_dim, temperature
export threading

const    outf="data.h5"
const    tebd_cutoff=1e-14
const    tebd_dt=0.05
const    tebd_order = 2  #2 or 4(4 doesn't seem to work as well as it should)
const    log_dt=0.1
const    final_time = 1.0
const    omega=1.0
const    t=1.0
const    gamma=sqrt(2.0)
const    threading="blocksparse"
const    N=4
const    boson_dim=21
const    temperature=0.4

end