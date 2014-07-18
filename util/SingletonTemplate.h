#ifndef __SINGLETON_HPP_
#define __SINGLETON_HPP_
#include <stddef.h>  // defines NULL
#include <cassert>  // defines NULL

template <class T>
class SingletonTemplate
{
public:
  static T& Instance() {
      if(!m_pInstance) m_pInstance = new T;
      assert(m_pInstance != NULL);
      return *m_pInstance;
  }
protected:
  SingletonTemplate();
  ~SingletonTemplate();
private:
  SingletonTemplate(SingletonTemplate const&);
  SingletonTemplate& operator=(SingletonTemplate const&);
  static T* m_pInstance;
};

template <class T> T* SingletonTemplate<T>::m_pInstance=NULL;

#endif
