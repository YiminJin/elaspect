/*
  Copyright (C) 2023 by Yimin Jin.

  This file is part of elASPECT.

  elASPECT is modified from the free software ASPECT; you can 
  redistribute it and/or modify it under the terms of the GNU 
  General Public License as published by the Free Software 
  Foundation; either version 2, or (at your option) any later 
  version.

  elASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with elASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _elaspect_plugins_h
#define _elaspect_plugins_h

#include <deal.II/base/utilities.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/exceptions.h>

#include <tuple>
#include <string>
#include <list>

namespace elaspect
{
  using namespace dealii;

  namespace Plugins
  {
    /**
     * This function returns if a given plugin (e.g. a material model returned
     * from SimulatorAccess::get_material_model() ) matches a certain plugin
     * type (e.g. MaterialModel::Simple). This check is needed, because often
     * it is only possible to get a reference to an Interface, not the actual
     * plugin type, but the actual plugin type might be important. For example
     * a radial gravity model might only be implemented for spherical geometry
     * models, and would want to check if the current geometry is in fact a
     * spherical shell.
     */
    template <typename TestType, typename PluginType>
    inline
    bool
    plugin_type_matches (const PluginType &object)
    {
      return (dynamic_cast<const TestType *> (&object) != nullptr);
    }

    /**
     * This function converts a reference to a type (in particular a reference
     * to an interface class) into a reference to a different type (in
     * particular a plugin class). This allows accessing members of the plugin
     * that are not specified in the interface class. Note that you should
     * first check if the plugin type is actually convertible by calling
     * plugin_matches_type() before calling this function. If the plugin is
     * not convertible this function throws an exception.
     */
    template <typename TestType, typename PluginType>
    inline
    TestType &
    get_plugin_as_type (PluginType &object)
    {
      AssertThrow(plugin_type_matches<TestType>(object),
                  ExcMessage("You have requested to convert a plugin of type <"
                             + boost::core::demangle(typeid(PluginType).name())
                             + "> into type <"
                             + boost::core::demangle(typeid(TestType).name()) +
                             ">, but this cast cannot be performed."));

      // We can safely dereference the pointer, because we checked above that
      // the object is actually of type TestType, and so the result
      // is not a nullptr.
      return *dynamic_cast<TestType *> (&object);
    }
  }

  namespace internal
  {
    namespace Plugins
    {
      template <typename InterfaceClass, 
                typename ModelClass>
      struct RegisterHelper
      {
        RegisterHelper (void (*register_function) (const std::string &,
                                                   const std::string &,
                                                   void ( *)(ParameterHandler &),
                                                   InterfaceClass * ( *)()),
                        const char *name,
                        const char *description)
        {
          register_function (name,
                             description,
                             &ModelClass::declare_parameters,
                             &factory);
        }

        static
        InterfaceClass *factory ()
        {
          return new ModelClass();
        }
      };

      template <typename InterfaceClass>
      struct PluginList
      {
        typedef
        std::tuple<std::string, 
                   std::string,
                   void ( *) (ParameterHandler &),
                   InterfaceClass *( *) ()> PluginInfo;

        static std::list<PluginInfo> *plugins;

        ~PluginList ();

        static
        void register_plugin (const std::string &name,
                              const std::string &description,
                              void (*declare_parameters_function) (ParameterHandler &),
                              InterfaceClass * (*factory_function) ());

        static
        std::string get_pattern_of_names ();

        static
        std::string get_description_string ();

        static
        void declare_parameters (ParameterHandler &prm);

        static
        InterfaceClass *
        create_plugin (const std::string &name,
                       const std::string &documentation);

        static
        InterfaceClass *
        create_plugin (const std::string  &name,
                       const std::string &documentation,
                       ParameterHandler &prm);

        DeclException1 (ExcUnknownPlugin,
                        std::string,
                        << "Can't create a plugin of name <" << arg1
                        << "> because such a plugin hasn't been declared.");
      };


      template <typename InterfaceClass>
      PluginList<InterfaceClass>::
      ~PluginList ()
      {
        // if any plugins have been registered, then delete
        // the list
        if (plugins != nullptr)
          delete plugins;
        plugins = nullptr;
      }


      template <typename InterfaceClass>
      void
      PluginList<InterfaceClass>::
      register_plugin (const std::string &name,
                       const std::string &description,
                       void (*declare_parameters_function) (ParameterHandler &),
                       InterfaceClass * (*factory_function) ())
      {
        // see if this is the first time we get into this
        // function and if so initialize the static member variable
        if (plugins == nullptr)
          plugins = new std::list<PluginInfo>();

        // verify that the same name has not previously been 
        // used to register a plugin, since we would then no
        // longer be able to identify the plugin
        for (typename std::list<PluginInfo>::const_iterator
             p = plugins->begin(); p != plugins->end(); ++p)
          Assert (std::get<0>(*p) != name,
                  ExcMessage ("A plugin with name <" + name + "> has "
                              "already been registered!"));

        // now add one record to the list
        plugins->push_back (PluginInfo(name,
                                       description,
                                       declare_parameters_function,
                                       factory_function));
      }


      template <typename InterfaceClass>
      std::string
      PluginList<InterfaceClass>::
      get_pattern_of_names ()
      {
        Assert (plugins != nullptr,
                ExcMessage ("No plugins registered!?"));

        // get all names and put them into a data structure that keeps
        // them sorted
        std::set<std::string> names;
        for (typename std::list<PluginInfo>::const_iterator
             p = plugins->begin(); p != plugins->end(); ++p)
          names.insert (std::get<0>(*p));

        // now create a pattern from all of these sorted names
        std::string pattern_of_names;
        for (typename std::set<std::string>::const_iterator
             p = names.begin();
             p != names.end(); ++p)
        {
          if (pattern_of_names.size() > 0)
            pattern_of_names += "|";
          pattern_of_names += *p;
        }

        return pattern_of_names;
      }


      template <typename InterfaceClass>
      std::string
      PluginList<InterfaceClass>::
      get_description_string ()
      {
        std::string description;

        // get all names_and_descriptions and put them into a data structure that keeps
        // them sorted
        std::map<std::string,std::string> names_and_descriptions;
        for (typename std::list<PluginInfo>::const_iterator
             p = plugins->begin(); p != plugins->end(); ++p)
          names_and_descriptions[std::get<0>(*p)] = std::get<1>(*p);;

        // then output it all
        typename std::map<std::string,std::string>::const_iterator
        p = names_and_descriptions.begin();
        while (true)
        {
          // write the name and
          // description of the
          // parameter
          description += "`";
          description += p->first;
          description += "': ";
          description += p->second;

          // increment the pointer
          // by one. if we are not
          // at the end yet then
          // add an empty line
          ++p;
          if (p != names_and_descriptions.end())
            description += "\n\n";
          else
            break;
        }

        return description;
      }


      template <typename InterfaceClass>
      void
      PluginList<InterfaceClass>::
      declare_parameters (ParameterHandler &prm)
      {
        Assert (plugins != nullptr,
                ExcMessage ("No postprocessors registered!?"));

        for (typename std::list<PluginInfo>::const_iterator
             p = plugins->begin();
             p != plugins->end(); ++p)
          (std::get<2>(*p))(prm);
      }


      template <typename InterfaceClass>
      InterfaceClass *
      PluginList<InterfaceClass>::
      create_plugin (const std::string &name,
                     const std::string &documentation)
      {
        (void)documentation;
        Assert (plugins != nullptr,
                ExcMessage ("No postprocessors registered!?"));
        AssertThrow (name != "unspecified",
                     ExcMessage(std::string("A plugin must have a name!\n\n"
                                            "This function was asked to create a plugin but no name for the "
                                            "plugin was provided. This may be due to the fact that you did not "
                                            "explicitly specify a name for this plugin in your input file and "
                                            "elASPECT does not provide a default for this kind of plugin, for "
                                            "example because no generally useful plugin exists. An example "
                                            "is that there is no default geometry: You need to explicitly "
                                            "provide one in the input file, and it seems like you have not "
                                            "done so.\n\n"
                                            "To find out which kind of plugin this function tries to create, "
                                            "take a look at the backtrace of this error message.\n\n"
                                            "The place that called this function also provided as "
                                            "additional information this:\n\n"
                                            "   <")
                                + documentation + ">"));

        for (typename std::list<PluginInfo>::const_iterator p = plugins->begin();
             p != plugins->end(); ++p)
          if (std::get<0>(*p) == name)
            {
              InterfaceClass *i = std::get<3>(*p)();
              return i;
            }

        AssertThrow (false, ExcUnknownPlugin(name));
        return nullptr;
      }


      template <typename InterfaceClass>
      InterfaceClass *
      PluginList<InterfaceClass>::
      create_plugin (const std::string &name,
                     const std::string &documentation,
                     ParameterHandler  &prm)
      {
        InterfaceClass *i = create_plugin(name, documentation);
        i->parse_parameters (prm);
        return i;
      }
    }
  }
}

#endif
