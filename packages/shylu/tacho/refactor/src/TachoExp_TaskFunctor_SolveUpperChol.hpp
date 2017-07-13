#ifndef __TACHOEXP_TASKFUNCTOR_SOLVE_UPPER_CHOL_HPP__
#define __TACHOEXP_TASKFUNCTOR_SOLVE_UPPER_CHOL_HPP__

/// \file TachoExp_TaskFunctor_FactorizeChol.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

#include "TachoExp_CholSupernodes.hpp"
#include "TachoExp_CholSupernodes_Serial.hpp"

namespace Tacho {
  
  namespace Experimental {
    
    template<typename MatValueType, typename ExecSpace>
    struct TaskFunctor_SolveUpperChol {
    public:
      typedef ExecSpace exec_space;

      typedef Kokkos::TaskScheduler<exec_space> sched_type;
      typedef typename sched_type::member_type member_type;
      
      typedef Kokkos::MemoryPool<exec_space> memory_pool_type;
      
      typedef int value_type; // functor return type
      typedef Kokkos::Future<int,exec_space> future_type;

      typedef MatValueType mat_value_type; // matrix value type

      typedef SupernodeInfo<mat_value_type,exec_space> supernode_info_type;
      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
      
    private:
      sched_type _sched;
      memory_pool_type _bufpool;

      supernode_info_type _info;
      ordinal_type _sid;
      
    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_SolveUpperChol() = delete;
      
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_SolveUpperChol(const sched_type &sched,
                                 const memory_pool_type &bufpool,
                                 const supernode_info_type &info,
                                 const ordinal_type sid)                                     
        : _sched(sched),
          _bufpool(bufpool),
          _info(info),
          _sid(sid) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        if (get_team_rank(member) == 0) {
          
          if (_info.serial_thres_size > _info.max_decendant_supernode_size(_sid)) {
            const ordinal_type n = _info.max_decendant_schur_size(_sid), nrhs = _info.x.dimension_1();
            const size_type bufsize = n*nrhs*sizeof(mat_value_type);
            
            mat_value_type *buf = bufsize > 0 ? (mat_value_type*)_bufpool.allocate(bufsize) : NULL;
            TACHO_TEST_FOR_ABORT(buf == NULL && bufsize != 0, "bufmemory pool allocation fails");   
            
            CholSupernodes<Algo::Workflow::Serial>
              ::solve_upper_recursive_serial(_sched, member, _info, _sid, true, buf, bufsize);
            
            _bufpool.deallocate(buf, bufsize);
          } else {
            // children information
            const ordinal_type 
              ibeg = _info.stree_ptr(_sid), 
              iend = _info.stree_ptr(_sid+1),
              isize = iend - ibeg;
            
            ordinal_type pm, pn; _info.getSuperPanelSize(_sid, pm, pn);
            const ordinal_type n = pn - pm, nrhs = _info.x.dimension_1();
            const size_type bufsize = n*nrhs*sizeof(mat_value_type);
            
            mat_value_type *buf = bufsize > 0 ? (mat_value_type*)_bufpool.allocate(bufsize) : NULL;
            TACHO_TEST_FOR_ABORT(buf == NULL && bufsize != 0, "bufmemory pool allocation fails");   
            
            CholSupernodes<Algo::Workflow::Serial>
              ::solve_upper_recursive_serial(_sched, member, _info, _sid, false, buf, bufsize);
            
            if (bufsize > 0)
              _bufpool.deallocate(buf, bufsize);
            
            // allocate dependence array to handle variable number of children schur contributions
            future_type depbuf[MaxDependenceSize] /* 3 */, *dep = &depbuf[0];
            const size_type depsize = (isize > MaxDependenceSize ? isize*sizeof(future_type) : 0);
            if (depsize) {
              dep = (future_type*)_sched.memory()->allocate(depsize);
              TACHO_TEST_FOR_ABORT(dep == NULL, "sched memory pool allocation fails");
              clear((char*)dep, depsize);
            }
            
            // spawn children tasks and this (their parent) depends on the children tasks
            for (ordinal_type i=0;i<isize;++i) {
              // the first child has a higher priority
              const auto priority = (i ? Kokkos::TaskPriority::Low : Kokkos::TaskPriority::High);
              
              const ordinal_type child = _info.stree_children(i+ibeg);
              auto f = Kokkos::task_spawn(Kokkos::TaskSingle(_sched, priority),
                                          TaskFunctor_SolveUpperChol(_sched, _bufpool, _info, child));
              TACHO_TEST_FOR_ABORT(f.is_null(), "task allocation fails");
              dep[i] = f;
            }
            
            // deallocate dependence array
            if (depsize) {
              // manually reset future to decrease the reference count
              for (ordinal_type i=0;i<isize;++i) (dep+i)->~future_type();
              _sched.memory()->deallocate((void*)dep, depsize);
            }
          }
        }
      }
    };
  }
}

#endif

            
            
