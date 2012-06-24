#ifndef UTIL_H
#define UTIL_H

#ifdef __GFORTRAN__
#define QQ(expr) "expr"
#else
#define QQ(expr) #expr
#endif

#define ASSERT(cond,msg) \
    if (.not. cond) call assert_msg(cond,QQ(cond),__FILE__,__LINE__,msg)

#define LOGMSG(verbosity,msg) \
  if (logv > verbosity .and. logv > 0) \
    write (F_LOG,'(2a,i0,2a)') __FILE__, ", line ", __LINE__, ": ", msg

#define WARN(str) \
  write(0,'(3a,i0,2a)') 'Warning: ', __FILE__, ", line ", __LINE__, ": ", str

#define WARN_IF(cond,msg,stmt) \
    if (.not. (cond)) then;; else; WARN(msg); stmt; endif

#endif
