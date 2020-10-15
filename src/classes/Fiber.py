import fibermodes

class fiber(object):

    def __init__(self,
                 core_radius,
                 core_index,
                 clad_radius,
                 clad_index):

        self.MaxDirect = 2 * clad_radius

        factory = fibermodes.FiberFactory()

        factory.addLayer(name='core',
                radius=core_radius,
                material='Fixed',
                geometry = "StepIndex",
                index=1.4489)

        factory.addLayer(name='cladding',
                         material='Fixed',
                         index=1)

        self.source = factory[0]
