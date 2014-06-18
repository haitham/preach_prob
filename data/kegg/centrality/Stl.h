//
//  Copyright 2012 Alin Dobra and Christopher Jermaine
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
#ifndef _STL_H_
#define _STL_H_

/** his file contains some macros to make life easier with stl
		containers and maps
*/


// Helper macros for Frequently Done Actions
// Macro to streamline scan of STL container
// list "list" is scanned and element is local variable that is
// instantiated as each element of the list
// Usage Scenario: 
//   FOREACH_STL(val, myList){
//     do something with val
//   }END_FOREACH
//
//
#define FOREACH_STL(el, list)																		\
	for(typeof(list.begin()) it = list.begin(); it != list.end(); it++){	\
	typeof(*it)& el = *it;
#ifndef END_FOREACH									
#define END_FOREACH }
#endif




#endif //  _STL_H_
