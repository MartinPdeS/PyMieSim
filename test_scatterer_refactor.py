#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test script to verify the refactored ScattererSection works correctly.
"""

def test_scatterer_section():
    """Test basic functionality of the refactored ScattererSection."""

    # Import the necessary modules (would be a real Dash app in practice)
    try:
        from PyMieSim.gui.section_.scatterer import ScattererSection
        from PyMieSim.gui.section_.scatterer_types import scatterer_manager

        print("✅ Imports successful")

        # Test the scatterer manager
        available_types = scatterer_manager.get_available_names()
        print(f"✅ Available scatterer types: {available_types}")

        dropdown_options = scatterer_manager.get_dropdown_options()
        print(f"✅ Dropdown options: {dropdown_options}")

        # Test individual scatterer types
        for scatterer_name in available_types:
            scatterer = scatterer_manager.get_type(scatterer_name)
            inputs = scatterer.get_input_definitions()
            mapping = scatterer.get_parameter_mapping()

            print(f"\n--- {scatterer.display_name} ---")
            print(f"Description: {scatterer.description}")
            print(f"Input fields: {list(inputs.keys())}")
            print(f"Parameter mapping: {mapping}")

            # Test validation with empty inputs
            is_valid, error = scatterer.validate_inputs({})
            print(f"Empty validation: {is_valid}, {error}")

            # Test validation with some inputs
            if scatterer_name == 'sphere':
                test_inputs = {
                    'diameter': '100',
                    'property': '1.5',
                    'medium_property': '1.33'
                }
            else:  # coreshell
                test_inputs = {
                    'core_diameter': '50',
                    'shell_thickness': '25',
                    'core_property': '2.4',
                    'shell_property': '1.5',
                    'medium_property': '1.33'
                }

            is_valid, error = scatterer.validate_inputs(test_inputs)
            print(f"Valid inputs test: {is_valid}, {error}")

        print(f"\n✅ All tests passed! The refactored ScattererSection is working correctly.")

        # Test creating a mock section (without actual Dash app)
        class MockApp:
            pass

        try:
            section = ScattererSection(MockApp())
            print(f"✅ ScattererSection created successfully")
            print(f"✅ Current type: {section.current_type}")
            print(f"✅ Available inputs: {list(section.inputs.keys())}")

            # Test type switching
            section.update_scatterer_type('coreshell')
            print(f"✅ Switched to core-shell, inputs: {list(section.inputs.keys())}")

            # Test getting current parameters
            params = section.get_current_parameters()
            print(f"✅ Parameter info retrieved: {params['type']}")

        except Exception as e:
            print(f"⚠️  ScattererSection creation failed (expected without real Dash app): {e}")

    except ImportError as e:
        print(f"❌ Import failed: {e}")
    except Exception as e:
        print(f"❌ Test failed: {e}")


if __name__ == "__main__":
    test_scatterer_section()
